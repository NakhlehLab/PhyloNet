package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.*;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class SymmetricDifference extends CommandBaseFileOut {

    private static final String _unnamedTaxonText = "[unnamed]";

    private NetworkNonEmpty _modelNetwork;

    private NetworkNonEmpty _experimentalNetwork;

    private boolean _contextChecked = false;

    public SymmetricDifference(SyntaxCommand directive, ArrayList<Parameter> params)
    {
        super(directive, params);
    }

    public boolean checkParams(Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String,Integer,Integer> errorDetected) {

         boolean noError = this.assertParamsCount(2, 3, errorDetected);

        if(noError)
        {
            return checkContext(sourceIdentToNetwork, errorDetected);
        }
        else
        {
            return noError;
        }
    }

    private boolean checkContext(Map<String,NetworkNonEmpty> sourceIdentToNetwork, Proc3<String,Integer,Integer> errorDetected)
    {
        boolean noError = true;

        Parameter modelTreeParam = params.get(0);
        noError = assertNetworkExists(sourceIdentToNetwork, modelTreeParam, errorDetected);

        String modelTreeParamValue = null;
        if(noError)
        {
            modelTreeParamValue = modelTreeParam.execute(GetSimpleParamValue.Singleton, null);
            _modelNetwork  = sourceIdentToNetwork.get(modelTreeParamValue);
        }

        Parameter experimentalTreeParam = params.get(1);
        noError = assertNetworkExists(sourceIdentToNetwork, experimentalTreeParam, errorDetected);

        String experimentalTreeParamValue = null;
        if(noError)
        {
            experimentalTreeParamValue = experimentalTreeParam.execute(GetSimpleParamValue.Singleton, null);
            _experimentalNetwork = sourceIdentToNetwork.get(experimentalTreeParamValue);
        }

        if(params.size() == 3)
        {
            Parameter outFileParam = params.get(2);
            noError = this.checkOutFileContext(outFileParam, errorDetected);
        }

        HashSet<Object> nonRootDegree1Model = _modelNetwork != null ? collectNonRootDegree1(_modelNetwork) : null;
        HashSet<Object> nonRootDegree1Exp   = _experimentalNetwork != null ? collectNonRootDegree1(_experimentalNetwork) : null;

       if(nonRootDegree1Model != null && nonRootDegree1Exp != null)
       {
           for(Object taxon : nonRootDegree1Model)
           {
               if(!nonRootDegree1Exp.remove(taxon))
               {
                   errorDetected.execute(String.format("Taxon '%s' in network '%s' does not appear in the network '%s'.",
                           taxon, modelTreeParamValue, experimentalTreeParamValue), modelTreeParam.getLine(), modelTreeParam.getColumn());
                   noError = false;
               }
           }

           for(Object taxon : nonRootDegree1Exp)
           {
               errorDetected.execute(String.format("Taxon '%s' in network '%s' does not appear in the network '%s'.",
                       taxon, experimentalTreeParamValue, modelTreeParamValue), experimentalTreeParam.getLine(), experimentalTreeParam.getColumn());
               noError = false;
           }
       }

        _contextChecked = true;
        return noError;
    }

    private HashSet<Object> collectNonRootDegree1(final NetworkNonEmpty network) {

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

    protected String produceResult() {
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

        return result;
    }


}
