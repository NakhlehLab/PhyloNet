package edu.rice.bioinfo.library.language.richnewick._1_0.csa;

import edu.rice.bioinfo.library.programming.Func;
import edu.rice.bioinfo.library.programming.Func1;

import javax.naming.NameNotFoundException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.LinkedList;

public class ContextAnalyser
{
    public static <SN,NN,E> CSAError[] Analyse(Iterable<SN> syntaxNodes,  final SyntaxNetworkInspector<SN> syntaxInspector,
                                               Iterable<NN> networkNodes, final NetworkInspector<NN,E> networkInspector,
                                               final Func1<NN,SN> getFirstSyntaxContributor)
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
             * Check to ensure bootstrap is a number between 0 and 1 inclusive.
             */
            String bootstrapStr = syntaxInspector.getBootstrapText(node);
            if(bootstrapStr != null)
            {
                try
                {
                    BigDecimal bootstrap = new BigDecimal(bootstrapStr, MathContext.UNLIMITED);

                    if(bootstrap.compareTo(BigDecimal.ZERO) == -1  || (bootstrap.compareTo(BigDecimal.ONE) == 1))
                    {
                        int bootstrapLine = syntaxInspector.getBootstrapLineNumber(node);
                        int bootstrapColumn = syntaxInspector.getBootstrapColumnNumber(node);

                        errors.addLast(
                                new CSAError(String.format("Bootstrap must be between zero and one inclusive.  Found '%s'.",
                                        bootstrapStr),bootstrapLine,bootstrapColumn));
                    }
                }
                catch(NumberFormatException e)
                {
                    int bootstrapLine = syntaxInspector.getBootstrapLineNumber(node);
                    int bootstrapColumn = syntaxInspector.getBootstrapColumnNumber(node);

                    errors.addLast(
                            new CSAError(String.format("Bootstrap '%s' is not a recognisable number.",
                                    bootstrapStr),bootstrapLine,bootstrapColumn));
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
                String parentProbText = networkInspector.getProbabilityText(edge);
                if(parentProbText != null)
                {
                    numParentsWithProb++;
                    BigDecimal parentProb = new BigDecimal(parentProbText, MathContext.UNLIMITED);
                    probSum = probSum.add(parentProb);
                }

                numParents++;
            }

            if(numParentsWithProb != numParents)
            {
                if(numParents != 1)
                {
                    SN syntaxNode = getFirstSyntaxContributor.execute(node);
                    String nodeLabel = syntaxInspector.getNodeLabelText(syntaxNode);

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
                }
            }
            else if(numParentsWithProb > 0 && probSum.compareTo(new BigDecimal("1", MathContext.UNLIMITED)) != 0)
            {
               SN syntaxNode = getFirstSyntaxContributor.execute(node);
               String nodeLabel = syntaxInspector.getNodeLabelText(syntaxNode);

                if(nodeLabel == null)
                {
                    errors.add(new CSAError("Unnamed node's parents' probabilities do not all sum to one but "  + probSum.toPlainString(), -1, -1));
                }
                else
                {
                   errors.add(new CSAError(String.format("Node %s's parents' probabilities do not all sum to one but "  + probSum.toPlainString(), nodeLabel),
                                           syntaxInspector.getNodeLabelTextLineNumber(syntaxNode),
                                           syntaxInspector.getNodeLabelTextColumnNumber(syntaxNode)));
                }
            }

        }

        return errors.toArray(new CSAError[0]);
    }
}
