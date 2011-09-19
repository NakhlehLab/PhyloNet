package edu.rice.cs.bioinfo.library.language.richnewick._1_0.csa;

import edu.rice.cs.bioinfo.library.programming.Func1;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

class ContextAnalyser
{
    public static <SN,NN,E> List<CSAError> analyse(Iterable<SN> syntaxNodes, final SyntaxNetworkInspector<SN> syntaxInspector,
                                                   Iterable<NN> networkNodes, final NetworkInspector<NN, E> networkInspector,
                                                   final Func1<NN, SN> getPrimarySyntaxContributor)
    {


        LinkedList<CSAError> errors = new LinkedList<CSAError>();
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
             * Check to ensure support is a number between 0 and 1 inclusive.
             */
            String supportStr = syntaxInspector.getSupportText(node);
            if(supportStr != null)
            {
                try
                {
                    BigDecimal support = new BigDecimal(supportStr, MathContext.UNLIMITED);

                    if(support.compareTo(BigDecimal.ZERO) == -1  || (support.compareTo(BigDecimal.ONE) == 1))
                    {
                        int supportLine = syntaxInspector.getSupportLineNumber(node);
                        int supportColumn = syntaxInspector.getSupportColumnNumber(node);

                        errors.addLast(
                                new CSAError(String.format("Support must be between zero and one inclusive.  Found '%s'.",
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
            * Check to ensure that, for all the in edges to a node, the probabilities sum to one.
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
                // probability assumed to be equal distirbution among all parents
                continue;
            }

            if(numParents == 1 && numParentsWithProb == 1 && probSum.compareTo(new BigDecimal("1", MathContext.UNLIMITED)) != 0)
            {
                if(nodeLabel == null)
                {
                    errors.add(new CSAError("Unnamed node with single parent must have probability 1.  Found  "  + probSum.toPlainString(), -1, -1));
                }
                else
                {
                    errors.add(new CSAError(String.format("Node %s with single parent must have probability 1.  Found "  + probSum.toPlainString(), nodeLabel),
                            syntaxInspector.getNodeLabelTextLineNumber(syntaxNode),
                            syntaxInspector.getNodeLabelTextColumnNumber(syntaxNode)));
                }

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

            if(probSum.compareTo(new BigDecimal("1", MathContext.UNLIMITED)) != 0)
            {
                if(nodeLabel == null)
                {
                    errors.add(new CSAError("Unnamed node's parents' probabilities do not all sum to 1. Found "  + probSum.toPlainString(), -1, -1));
                }
                else
                {
                    errors.add(new CSAError(String.format("Node %s's parents' probabilities do not all sum to 1.  Found "  + probSum.toPlainString(), nodeLabel),
                            syntaxInspector.getNodeLabelTextLineNumber(syntaxNode),
                            syntaxInspector.getNodeLabelTextColumnNumber(syntaxNode)));
                }
            }

        }

        return Collections.unmodifiableList(errors);
    }
}
