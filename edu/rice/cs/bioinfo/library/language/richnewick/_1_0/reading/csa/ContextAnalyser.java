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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.csa;

import edu.rice.cs.bioinfo.library.language.richnewick.reading.csa.CSAError;
import edu.rice.cs.bioinfo.library.programming.Func1;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

public class ContextAnalyser
{
    private final BigDecimal _hybridSumTolerance;

    public ContextAnalyser(BigDecimal hybridSumTolerance)
    {
        _hybridSumTolerance = hybridSumTolerance;
    }

    public ContextAnalyser()
    {
        this(BigDecimal.ZERO);
    }

    public <SN,NN,E> List<CSAError> analyse(Iterable<SN> syntaxNodes, final SyntaxNetworkInspector<SN> syntaxInspector,
                                                   Iterable<NN> networkNodes, final NetworkInspector<NN, E> networkInspector,
                                                   final Func1<NN, SN> getPrimarySyntaxContributor, boolean isRooted)
    {

        LinkedList<CSAError> errors = new LinkedList<CSAError>();
        if(!isRooted)
        {
            NN rootNode = networkInspector.getRootNode();
            SN rootNodeSyntax = getPrimarySyntaxContributor.execute(rootNode);
            String rootLabel = syntaxInspector.getNodeLabelText(rootNodeSyntax);

            int rootOutEdgesCount = 0;

            for(E edge : networkInspector.getEdges())
            {
                if(networkInspector.getTail(edge).equals(rootNode))
                {
                    rootOutEdgesCount++;
                }
            }

            if(rootLabel != null && rootOutEdgesCount == 2)
            {
                errors.add(new CSAError("Unrooted trees of the form ((...)A, (...)B) may not contain labes at the root position because A and B are considered adjacent.",
                                         syntaxInspector.getNodeLabelTextLineNumber(rootNodeSyntax), syntaxInspector.getNodeLabelTextColumnNumber(rootNodeSyntax)));
            }
        }


        for(SN node : syntaxNodes)
        {
            /*
             * Check to ensure branch length is a number.
             */
            String branchLengthStr = syntaxInspector.getBranchLengthText(node);
            if(branchLengthStr != null)
            {
                try
                {
                    Double branchLength = new Double(branchLengthStr);
                }
                catch(NumberFormatException e)
                {
                    int branchLengthLine = syntaxInspector.getBranchLengthLineNumber(node);
                    int branchLengthColumn = syntaxInspector.getBranchLengthColumnNumber(node);

                    errors.addLast(
                            new CSAError(String.format("Branch length '%s' is not a recognisable number.",
                                    branchLengthStr),branchLengthLine,branchLengthColumn));
                }
            }

            /*
             * Check to ensure support is a number between 0 and 100 inclusive.
             */
            String supportStr = syntaxInspector.getSupportText(node);
            if(supportStr != null)
            {
                try
                {
                    BigDecimal support = new BigDecimal(supportStr, MathContext.UNLIMITED);

                    if(support.compareTo(BigDecimal.ZERO) == -1  || (support.compareTo(BigDecimal.valueOf(100)) == 1))
                    {
                        int supportLine = syntaxInspector.getSupportLineNumber(node);
                        int supportColumn = syntaxInspector.getSupportColumnNumber(node);

                        errors.addLast(
                                new CSAError(String.format("Support must be between zero and a hundred inclusive.  Found '%s'.",
                                        supportStr),supportLine,supportColumn));
                    }
                }
                catch(NumberFormatException e)
                {
                    int supportLine = syntaxInspector.getSupportLineNumber(node);
                    int supportColumn = syntaxInspector.getSupportColumnNumber(node);

                    errors.addLast(
                            new CSAError(String.format("Support '%s' is not a recognisable number.",
                                    supportStr),supportLine,supportColumn));
                }
            }

            /*
             * Check to ensure probability is a number between 0 and 1 inclusive.
             */
            String probStr = syntaxInspector.getProbabilityText(node);
            if(probStr != null)
            {
                try
                {
                    BigDecimal prob = new BigDecimal(probStr, MathContext.UNLIMITED);

                    if(prob.compareTo(BigDecimal.ZERO) == -1  || (prob.compareTo(BigDecimal.ONE) == 1))
                    {
                        int probLine = syntaxInspector.getProbabilityLineNumber(node);
                        int probColumn = syntaxInspector.getProbabilityColumnNumber(node);

                        errors.addLast(
                                new CSAError(String.format("Probability must be between zero and one inclusive.  Found '%s'.",
                                        probStr),probLine,probColumn));
                    }
                }
                catch(NumberFormatException e)
                {
                    int probLine = syntaxInspector.getProbabilityLineNumber(node);
                    int probColumn = syntaxInspector.getProbabilityColumnNumber(node);

                    errors.addLast(
                            new CSAError(String.format("Probability '%s' is not a recognisable number.",
                                    probStr),probLine,probColumn));
                }
            }

            /*
             * Check to ensure hybrid node type is R, H or LGT
             */
            String hybridNodeType = syntaxInspector.getHybridNodeType(node);
            if(hybridNodeType != null)
            {
                if(!hybridNodeType.equals("R") && ! hybridNodeType.equals("H") &&
                        !hybridNodeType.equals("LGT"))
                {

                    int hybridLine = syntaxInspector.getHybridNodeTypeLineNumber(node);
                    int hybridColumn = syntaxInspector.getHybridNodeTypeColumnNumber(node);

                    errors.addLast(
                            new CSAError(String.format("Hybrid node type '%s' is unknown. Expected 'R', 'H', or 'LGT'.",
                                    hybridNodeType),hybridLine,hybridColumn));

                }
            }

            /*
             * Check to ensure hybrid node index is greater than or equal to one
             */
            String hybridNodeIndex = syntaxInspector.getHybridNodeIndexText(node);
            if(hybridNodeIndex != null)
            {
                try
                {
                    Integer index = new Integer(hybridNodeIndex);

                    if(index < 1)
                    {
                        int indexLine = syntaxInspector.getHybridNodeIndexLineNumber(node);
                        int indexColumn = syntaxInspector.getHybridNodeIndexColumnNumber(node);

                        errors.addLast(
                                new CSAError(String.format("Hybrid node index must be greater than or equal to one.  Found '%s'.",
                                        hybridNodeIndex),indexLine,indexColumn));
                    }
                }
                catch(NumberFormatException e)
                {
                    int indexLine = syntaxInspector.getHybridNodeIndexLineNumber(node);
                    int indexColumn = syntaxInspector.getHybridNodeIndexColumnNumber(node);

                    errors.addLast(
                            new CSAError(String.format("Hybrid node index '%s' is not a recognisable number.",
                                    hybridNodeIndex),indexLine,indexColumn));
                }
            }


        }

        for(NN node : networkNodes)
        {
            /*
            * Check to ensure that, for all the in edges to a node, the probabilities sum to one within tolerance.
            */

            BigDecimal probSum = new BigDecimal("0", MathContext.UNLIMITED);
            int numParentsWithProb = 0;
            int numParents = 0;
            for(E edge : networkInspector.getAllInEdges(node))
            {
                String parentProbText = networkInspector.getEdgeProbabilityText(edge);
                if(parentProbText != null)
                {
                    numParentsWithProb++;
                    BigDecimal parentProb = new BigDecimal(parentProbText, MathContext.UNLIMITED);
                    probSum = probSum.add(parentProb);
                }

                numParents++;
            }

            SN syntaxNode = getPrimarySyntaxContributor.execute(node);
            String nodeLabel = syntaxInspector.getNodeLabelText(syntaxNode);

            if(numParentsWithProb == 0)
            {
                // probability assumed to be equal distribution among all parents
                continue;
            }

            if(numParentsWithProb != numParents)
            {


                if(nodeLabel == null)
                {
                    errors.add(new CSAError("Unnamed node's parents' do not all contain probabilities.", -1, -1));
                }
                else
                {
                    errors.add(new CSAError(String.format("Node %s's parents' do not all contain probabilities.", nodeLabel),
                            syntaxInspector.getNodeLabelTextLineNumber(syntaxNode),
                            syntaxInspector.getNodeLabelTextColumnNumber(syntaxNode)));
                }
                continue;
            }

            if(probSum.subtract(BigDecimal.ONE).abs().compareTo(_hybridSumTolerance) == 1)
            {
                if(nodeLabel == null)
                {
                    errors.add(new CSAError("Unnamed node's parents' probabilities do not all sum to 1 within a tolerance of " + _hybridSumTolerance + ". Found "  + probSum.toPlainString(), -1, -1));
                }
                else
                {
                    errors.add(new CSAError(String.format("Node %s's parents' probabilities do not all sum to 1 within a tolerance of " + _hybridSumTolerance + ". Found " + probSum.toPlainString(), nodeLabel),
                            syntaxInspector.getNodeLabelTextLineNumber(syntaxNode),
                            syntaxInspector.getNodeLabelTextColumnNumber(syntaxNode)));
                }
            }

        }

        return Collections.unmodifiableList(errors);
    }
}
