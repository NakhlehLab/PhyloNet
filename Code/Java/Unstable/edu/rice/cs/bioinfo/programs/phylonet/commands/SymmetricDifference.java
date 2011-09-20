package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.*;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.awt.geom.Path2D;
import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class SymmetricDifference extends CommandBase {

    private static final String _unnamedTaxonText = "[unnamed]";

    private final ArrayList<Parameter> _params;

    private Network _modelNetwork;

    private Network _experimentalNetwork;

    private boolean _contextChecked = false;

    private File _outFile;

    public SymmetricDifference(SyntaxCommand directive, ArrayList<Parameter> params)
    {
        super(directive);
        _params = params;
    }

    public boolean checkParams(Proc3<String,Integer,Integer> errorDetected) {

          boolean noError = true;

          if(2 > _params.size())
           {
               SyntaxCommand directive = this.getDefiningSyntaxCommand();
                errorDetected.execute(String.format("Expected at least 2 parameters for command %s but found %s.",
                                                    directive.getName(), _params.size()),
                                                    directive.getLine(), directive.getColumn());
               noError = false;
           }

           if(3 < _params.size())
           {
                SyntaxCommand directive = this.getDefiningSyntaxCommand();
                errorDetected.execute(String.format("Expected at most 3 parameters for command %s but found %s.",
                                                    directive.getName(), _params.size()),
                                                    directive.getLine(), directive.getColumn());
               noError = false;
           }

        return noError;
    }

    public boolean checkContext(Map<String,Network> sourceIdentToNetwork, Proc3<String,Integer,Integer> errorDetected)
    {
         boolean noError = true;

        Parameter modelTreeParam = _params.get(0);
        if(sourceIdentToNetwork.containsKey(modelTreeParam.getValue()))
        {
            _modelNetwork  = sourceIdentToNetwork.get(modelTreeParam.getValue());
        }
        else
        {
            noError = false;
            errorDetected.execute(String.format("Unknown identifier '%s'.", modelTreeParam.getValue()),
                                                modelTreeParam.getLine(), modelTreeParam.getColumn());
        }

        Parameter experimentalTreeParam = _params.get(1);
        if(sourceIdentToNetwork.containsKey(experimentalTreeParam.getValue()))
        {
            _experimentalNetwork = sourceIdentToNetwork.get(experimentalTreeParam.getValue());
        }
        else
        {
            noError = false;
            errorDetected.execute(String.format("Unknown identifier '%s'.", experimentalTreeParam.getValue()),
                                                experimentalTreeParam.getLine(), experimentalTreeParam.getColumn());
        }

        if(_params.size() == 3)
        {
            Parameter outFileParam = _params.get(2);
            _outFile = new File(outFileParam.getValue());

            if(!_outFile.exists())
            {
                try
                {
                    _outFile.createNewFile();
                    _outFile.delete();
                }
                catch(IOException e)
                {
                     noError = false;
                    errorDetected.execute(String.format("Invalid file name: '%s'.", outFileParam.getValue()),
                                                         outFileParam.getLine(), outFileParam.getColumn());
                }
                catch(SecurityException e)
                {
                     noError = false;
                     errorDetected.execute(String.format("No access to file: '%s'.", outFileParam.getValue()),
                                                         outFileParam.getLine(), outFileParam.getColumn());
                }
            }
        }

        HashSet<Object> nonRootDegree1Model = collectNonRootDegree1(_modelNetwork);
        HashSet<Object> nonRootDegree1Exp = collectNonRootDegree1(_experimentalNetwork);

        NodeLabelAlgo<String,Object,RuntimeException> taxonNameExtractor = new NodeLabelAlgo<String, Object, RuntimeException>() {
                    public String forNodeLabelNonEmpty(NodeLabelNonEmpty nodeLabelNonEmpty, Object o) throws RuntimeException {
                        return nodeLabelNonEmpty.Label.Content;
                    }

                    public String forNodeLabelEmpty(NodeLabelEmpty nodeLabelEmpty, Object o) throws RuntimeException {
                        return _unnamedTaxonText;
                    }
                };

        for(Object taxon : nonRootDegree1Model)
        {
            if(!nonRootDegree1Exp.remove(taxon))
            {
                errorDetected.execute(String.format("Taxon '%s' in network '%s' does not appear in the network '%s'. ",
                                                    taxon, modelTreeParam.getValue(), experimentalTreeParam.getValue()), modelTreeParam.getLine(), modelTreeParam.getColumn());
                noError = false;
            }
        }

        for(Object taxon : nonRootDegree1Exp)
        {
            errorDetected.execute(String.format("Taxon '%s' in network '%s' does not appear in the network '%s'. ",
                                                taxon, experimentalTreeParam.getValue(), modelTreeParam.getValue()), experimentalTreeParam.getLine(), experimentalTreeParam.getColumn());
            noError = false;
        }

        _contextChecked = true;
        return noError;
    }

    private HashSet<Object> collectNonRootDegree1(final Network network) {

        final HashSet<Object> tbr = new HashSet<Object>();

        int numSubtreesOfPrinciple = 0;
        for(Subtree tree : network.PrincipleDescendants.Subtrees)
        {
            numSubtreesOfPrinciple++;
            collectNonRootDegree1Help(tree.NetworkInfo, tree.Descendants, tbr);
        }

        if(numSubtreesOfPrinciple == 1)
        {
            network.RootageQualifier.execute(new RootageQualifierAlgo<Object, Object, RuntimeException>() {
                public Object forEmptyQualifier(RootageQualifierEmpty rootageQualifierEmpty, Object o) throws RuntimeException {
                    return null;
                }

                public Object forNonEmptyQualifier(RootageQualifierNonEmpty rootageQualifierNonEmpty, Object o) throws RuntimeException {

                    if(!rootageQualifierNonEmpty.isRooted())
                    {
                        tbr.add(makeTaxonObject(network.PrincipleInfo));
                    }

                    return null;

                }
            }, null);
        }

        return tbr;


    }



    private void collectNonRootDegree1Help(NetworkInfo parent, DescendantList descendents, HashSet<Object> nodesAccum)
    {
        int numSubtrees = 0;
        for(Subtree tree : descendents.Subtrees)
        {
            numSubtrees++;
            collectNonRootDegree1Help(tree.NetworkInfo, tree.Descendants, nodesAccum);
        }

        if(numSubtrees == 0)
        {
            nodesAccum.add(makeTaxonObject(parent));
        }
    }

    private Object makeTaxonObject(NetworkInfo node)
    {
        return node.NodeLabel.execute(new NodeLabelAlgo<Object, Object, RuntimeException>() {
                public Object forNodeLabelNonEmpty(NodeLabelNonEmpty nodeLabelNonEmpty, Object o) throws RuntimeException {

                    return nodeLabelNonEmpty.Label.Content;
                }

                public Object forNodeLabelEmpty(NodeLabelEmpty nodeLabelEmpty, Object o) throws RuntimeException {

                   return new Object()
                   {
                       @Override public String toString()
                       {
                           return _unnamedTaxonText;
                       }
                   };
                }
            }, null);
    }

    public void executeCommand(Proc<String> displayResult) throws IOException {
        if(!_contextChecked)
        {
            throw new IllegalStateException("checkContext must be called prior to execute.");
        }

        Tree modelTree = NetworkTransformer.toTree(_modelNetwork);
        Tree experamentalTree = NetworkTransformer.toTree(_experimentalNetwork);

        String m = modelTree.toStringWD();
        String e = experamentalTree.toStringWD();

        final edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference sd =
                new edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference();

        sd.computeDifference(experamentalTree, modelTree);


        String result =
                String.format("False Negatives: %s\nFalse Positives: %s\n# Internal Edges Model: %s\n# Internal Edges Experimental: %s\n",
                        sd.getFalseNegativeCount(), sd.getFalsePositiveCount(), sd.getNumInternalEdges1(), sd.getNumInternalEdges2());

        if(_outFile == null)
        {
            displayResult.execute(result);
        }
        else
        {
            if(_outFile.exists())
            {
                _outFile.delete();

            }
            _outFile.createNewFile();
            FileWriter fw = new FileWriter(_outFile);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write(result);
            bw.flush();
            fw.flush();
            bw.close();
            fw.close();
        }
    }


    public <R, T, E extends Exception> R execute(CommandAlgo<R, T, E> algo, T input) throws E {
        return algo.forSymmetricDifference(this, input);
    }
}
