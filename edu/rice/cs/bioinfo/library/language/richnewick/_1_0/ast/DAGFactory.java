package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilder;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.HashMap;
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

            public Object forNetworkNonEmpty(NetworkNonEmpty network, Object input) throws RuntimeException {

                 Map<BigInteger,T> hybridIndexToNode = new HashMap<BigInteger, T>();
                 T parent = createNode(network.PrincipleInfo, builder, hybridIndexToNode);
                 processDescendantList(parent, network.PrincipleDescendants, builder, hybridIndexToNode);
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
