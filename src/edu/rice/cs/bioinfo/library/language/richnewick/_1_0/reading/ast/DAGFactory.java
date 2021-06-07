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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.GraphBuilder;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/26/11
 * Time: 12:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class DAGFactory
{
    public static <T> void makeDAG(Network network, final GraphBuilder<T> builder)
    {
        network.execute(new NetworkAlgo<Object, Object, RuntimeException>() {
            public Object forNetworkEmpty(NetworkEmpty network, Object input) throws RuntimeException {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public Object forNetworkNonEmpty(final NetworkNonEmpty network, Object input) throws RuntimeException {

                final Map<BigInteger,T> hybridIndexToNode = new HashMap<BigInteger, T>();

                network.RootageQualifier.execute(new RootageQualifierAlgo<Object, Object, RuntimeException>() {
                    public Object forEmptyQualifier(RootageQualifierEmpty rootage, Object input) throws RuntimeException {


                        performBuild(network.PrincipleInfo, network.PrincipleDescendants);
                        return null;
                    }

                    public Object forNonEmptyQualifier(RootageQualifierNonEmpty rootage, Object input) throws RuntimeException {

                        if(!rootage.isRooted())
                        {
                            // unrooted trees where root has only two children maps to special graph topology.
                            Iterator<Subtree> principleSubsIt = network.PrincipleDescendants.Subtrees.iterator();

                            if(principleSubsIt.hasNext())
                            {
                               Subtree syntaxRootFirstChild = principleSubsIt.next();
                               if(principleSubsIt.hasNext())
                               {
                                  Subtree syntaxRootSecondChild = principleSubsIt.next();

                                   if(!principleSubsIt.hasNext()) // unrooted and two children from syntax root, speical case
                                   {
                                       LinkedList<Subtree> newFirstChildChildren = new LinkedList<Subtree>();

                                       for(Subtree firstChildChild : syntaxRootFirstChild.Descendants.Subtrees)
                                       {
                                            newFirstChildChildren.add(firstChildChild);
                                       }
                                       newFirstChildChildren.add(syntaxRootSecondChild);

                                       DescendantList newDecList = new DescendantList(newFirstChildChildren);

                                       performBuild(syntaxRootFirstChild.NetworkInfo, newDecList);
                                       return null;

                                   }
                               }
                            }
                        }

                        performBuild(network.PrincipleInfo, network.PrincipleDescendants);
                        return null;

                    }

                    private void performBuild(NetworkInfo principleInfo, DescendantList pincipleDesc)
                    {
                        T parent = createNode(principleInfo, builder, hybridIndexToNode);
                        processDescendantList(parent, pincipleDesc, builder, hybridIndexToNode);

                    }

                }, null);








                return null;
            }
        }, null);


    }


    private static <T> void processDescendantList(T parent, DescendantList children, GraphBuilder<T> builder,
                                                  Map<BigInteger,T> hybridIndexToNode)
    {
        for(Subtree sub : children.Subtrees)
        {
            BigDecimal branchLength = sub.NetworkInfo.BranchLength.execute(new BranchLengthAlgo<BigDecimal, Object, RuntimeException>() {
                public BigDecimal forBranchLengthEmpty(BranchLengthEmpty branchLength, Object input) {
                    return null;
                }

                public BigDecimal forBranchLengthNonEmpty(BranchLengthNonEmpty branchLength, Object input)  {
                    return new BigDecimal(branchLength.Length.Content);
                }
            }, null);

            BigDecimal support = sub.NetworkInfo.Support.execute(new SupportAlgo<BigDecimal, Object, RuntimeException>() {
                public BigDecimal forSupportNonEmpty(SupportNonEmpty support, Object input)  {
                    return new BigDecimal(support.SupportValue.Content);
                }

                public BigDecimal forSupportEmpty(SupportEmpty support, Object input)  {
                    return null;
                }
            }, null);

            BigDecimal probability = sub.NetworkInfo.Probability.execute(new ProbabilityAlgo<BigDecimal, Object, RuntimeException>() {
                public BigDecimal forProbabilityEmpty(ProbabilityEmpty prob, Object input)  {
                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public BigDecimal forProbabilityNonEmpty(ProbabilityNonEmpty prob, Object input)  {
                    return new BigDecimal(prob.ProbabilityValue.Content);
                }
            }, null);

            T child = createNode(sub.NetworkInfo, builder, hybridIndexToNode);
            builder.createDirectedEdge(parent, child, branchLength, support, probability);
            processDescendantList(child, sub.Descendants, builder, hybridIndexToNode);
        }
    }

    private static <T> T createNode(final NetworkInfo nodeInfo, final GraphBuilder<T> builder,
                                    final Map<BigInteger,T> hybridIndexToNode)
    {
        final String nodeLabel = nodeInfo.NodeLabel.execute(new NodeLabelAlgo<String, Object, RuntimeException>() {

                    public String forNodeLabelNonEmpty(NodeLabelNonEmpty node, Object input) {
                        return node.Label.Content;
                    }

                    public String forNodeLabelEmpty(NodeLabelEmpty node, Object input) {
                        return null;
                    }
                }, null);

        return nodeInfo.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<T, Object, RuntimeException>() {

            public T forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input) {
                return builder.createNode(nodeLabel);
            }

            public T forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input) {

               return onNotEmpty(qualifier, HybridNodeType.Unspecified);
            }

            public T forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input) {
                return onNotEmpty(qualifier, HybridNodeType.fromString(qualifier.HybridNodeType.Content));
            }

            private T onNotEmpty(HybridNodeQualifierNonEmpty qualifier, HybridNodeType type)
            {
                 BigInteger hybridIndex = new BigInteger(qualifier.HybridNodeIndex.Content);

                if(hybridIndexToNode.containsKey(hybridIndex))
                {
                    return hybridIndexToNode.get(hybridIndex);
                }
                else
                {
                    T node = builder.createHybridNode(nodeLabel, type , hybridIndex);
                    hybridIndexToNode.put(hybridIndex, node);
                    return node;
                }
            }

        }, null);


    }

}
