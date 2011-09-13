package edu.rice.bioinfo.programs.phylonet.commands;

import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.Parameter;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SyntaxCommand;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.*;
import edu.rice.bioinfo.library.programming.Proc3;
import edu.rice.bioinfo.library.programming.extensions.java.lang.iterable.Proc;
import edu.rice.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.awt.font.ShapeGraphicAttribute;
import java.io.FileOutputStream;
import java.io.PrintStream;
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
public class SymmetricDifference implements Command {

    private static final String _unnamedTaxonText = "[unnamed]";

    public interface SymmetricDifferenceResult
    {
        public int getFalseNegativesCount();

        public int getFalsePositivesCount();

        public int getModelInternalEdgesCount();

        public int getExperamentalInternalEdgesCount();
    }

    public void checkParams(SyntaxCommand command, ArrayList<Parameter> params, Proc3<String,Integer,Integer> errorDetected) {


          if(2 > params.size())
           {
                errorDetected.execute(String.format("Expected at least 2 parameters for command %s but found %s.",
                                                    command.getName(), params.size()),
                                                    command.getLine(), command.getColumn());
           }

           if(3 < params.size())
           {
                errorDetected.execute(String.format("Expected at most 3 parameters for command %s but found %s.",
                                                    command.getName(), params.size()),
                                                    command.getLine(), command.getColumn());
           }
    }

    public boolean checkContext(Network modelNetwork, String modelNetworkIdent, Network experimentalNetwork, String experimetnalNetworkIdent, Proc<String> errorDetected)
    {
        boolean noError = true;

        HashSet<Object> nonRootDegree1Model = collectNonRootDegree1(modelNetwork);
        HashSet<Object> nonRootDegree1Exp = collectNonRootDegree1(experimentalNetwork);

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
                errorDetected.execute(String.format("Taxon '%s' in network '%s' does not appear in the network '%s'. ", taxon, modelNetworkIdent, experimetnalNetworkIdent));
                noError = false;
            }
        }

        for(Object taxon : nonRootDegree1Exp)
        {
            errorDetected.execute(String.format("Taxon '%s' in network '%s' does not appear in the network '%s'. ", taxon, experimetnalNetworkIdent, modelNetworkIdent));
            noError = false;
        }

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

    public SymmetricDifferenceResult execute(Network modelNetwork, Network experamentalNetwork)
    {
        Tree modelTree = NetworkTransformer.toTree(modelNetwork);
        Tree experamentalTree = NetworkTransformer.toTree(experamentalNetwork);

        String m = modelTree.toStringWD();
        String e = experamentalTree.toStringWD();

        final edu.rice.bioinfo.programs.phylonet.algos.SymmetricDifference sd =
                new edu.rice.bioinfo.programs.phylonet.algos.SymmetricDifference();

        sd.computeDifference(experamentalTree, modelTree);

        return new SymmetricDifferenceResult() {
            public int getFalseNegativesCount() {
                return sd.getFalseNegativeCount();
            }

            public int getFalsePositivesCount() {
                return sd.getFalsePositiveCount();
            }

            public int getModelInternalEdgesCount() {
               return sd.getNumInternalEdges1();
            }

            public int getExperamentalInternalEdgesCount() {
                return sd.getNumInternalEdges2();
            }
        };


    }


    public <R, T, E extends Exception> R execute(CommandAlgo<R, T, E> algo, T input) throws E {
        return algo.forSymmetricDifference(this, input);
    }
}
