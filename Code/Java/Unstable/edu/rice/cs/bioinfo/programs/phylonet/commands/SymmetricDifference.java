/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

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
@CommandName("SymmetricDifference")
public class SymmetricDifference extends CommandBaseFileOut {

    private static final String _unnamedTaxonText = "[unnamed]";

    private NetworkNonEmpty _modelNetwork;

    private NetworkNonEmpty _experimentalNetwork;

    private boolean _rooted = false;

    public SymmetricDifference(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                               Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    protected int getMinNumParams()
    {
        return 2;
    }

    protected int getMaxNumParams()
    {
        return 3;
    }

    public boolean checkParamsForCommand() {

        boolean noError = true;

        ParamExtractor rParam = new ParamExtractor("r", this.params, this.errorDetected);
        if(rParam.ContainsSwitch)
        {
            _rooted = true;
        }


        noError = noError && checkForUnknownSwitches("r");
        checkAndSetOutFile(rParam);

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
        _modelNetwork = this.assertAndGetTree(0);
        noError = noError && _modelNetwork != null;
        String modelTreeParamValue = _modelNetwork == null ? null : modelTreeParam.execute(GetSimpleParamValue.Singleton, null);


        Parameter experimentalTreeParam = params.get(1);
        _experimentalNetwork = this.assertAndGetTree(1);
        noError = noError && _experimentalNetwork != null;
        String experimentalTreeParamValue = _experimentalNetwork == null ? null : experimentalTreeParam.execute(GetSimpleParamValue.Singleton, null);


        if(params.size() == 4)
        {
            noError = this.checkOutFileContext(3);
        }
        if(params.size() == 3 && _rooted == false)
        {
            noError = this.checkOutFileContext(2);
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


        Tree modelTree = NetworkTransformer.toTree(_modelNetwork);
        Tree experimentalTree = NetworkTransformer.toTree(_experimentalNetwork);

        final edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference sd =
                new edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference();

        sd.computeDifference(experimentalTree, modelTree, _rooted);

        //String result = String.format("\nFalse Negatives: %s\nFalse Positives: %s\n# Internal Edges Model: %s\n# Internal Edges Experimental: %s",sd.getFalseNegativeCount(), sd.getFalsePositiveCount(), sd.getNumInternalEdges1(), sd.getNumInternalEdges2());
        String result =
                String.format("\n# False Positive Edges: %s\n# False Negative Edges: %s\n# Internal Edges Model: %s\n# Internal Edges Experimental: %s" +
                        "\nNormalized False Positive: %s\nNormalized False Negative: %s\nNormalized RF-Distance: %s",
                        sd.getFalsePositiveCount(), sd.getFalseNegativeCount(), sd.getNumInternalEdges2(), sd.getNumInternalEdges1(),
                        sd.getWeightedFalsePositive(), sd.getWeightedFalseNegative(), sd.getWeightedAverage());


        return result;
    }


}
