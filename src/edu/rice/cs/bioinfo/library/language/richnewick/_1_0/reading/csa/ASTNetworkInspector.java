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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/28/11
 * Time: 1:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class ASTNetworkInspector implements SyntaxNetworkInspector<NetworkInfo>, NetworkInspector<Object, ASTNetworkInspector.Edge>
{
    public class Edge
    {
        public final Object Tail;

        public final Object Tip;

        public final String ProbabilityText;

        public Edge(Object tail, Object tip, String probabilityText)
        {
            Tail = tail;
            Tip = tip;
            ProbabilityText = probabilityText;
        }
    }

    private Map<Object, LinkedList<Edge>> _networkNodeToInEdges = new HashMap<Object, LinkedList<Edge>>();

    private LinkedList<Edge> _edges = new LinkedList<Edge>();

    public Iterable<Edge> getEdges()
    {
        return _edges;
    }

    private LinkedList<NetworkInfo> _syntaxNodes = new LinkedList<NetworkInfo>();

    private Map<Object, NetworkInfo> _networkNodeToPrimarySyntaxNode = new HashMap<Object, NetworkInfo>();

    public Iterable<NetworkInfo> getSyntaxNodes()
    {
        return _syntaxNodes;
    }

    private LinkedList<Object> _networkNodes = new LinkedList<Object>();

    public Iterable<Object> getNetworkNodes()
    {
        return _networkNodes;
    }

    private final NetworkInfo _rootNode;

    public NetworkInfo getRootNode()
    {
        return _rootNode;
    }

    public ASTNetworkInspector(Network network)
    {
         _rootNode = network.execute(new NetworkAlgo<NetworkInfo, Object, RuntimeException>() {
            public NetworkInfo forNetworkEmpty(NetworkEmpty network, Object input) throws RuntimeException {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public NetworkInfo forNetworkNonEmpty(NetworkNonEmpty network, Object input) throws RuntimeException {
                gatherInEdges(network.PrincipleInfo, network.PrincipleDescendants);
                _networkNodes.add(network.PrincipleInfo);
                _networkNodeToInEdges.put(network.PrincipleInfo, new LinkedList<Edge>());
                _networkNodeToPrimarySyntaxNode.put(network.PrincipleInfo, network.PrincipleInfo);
                return network.PrincipleInfo;
            }
        }, null);


    }

    private void gatherInEdges(final NetworkInfo parent, DescendantList descendants)
    {
        _syntaxNodes.add(parent);

        for(Subtree subTree : descendants.Subtrees)
        {
            final NetworkInfo child =  subTree.NetworkInfo;

            child.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<Object, Object, RuntimeException>() {

                public Object forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input)  {

                    Edge edge = new Edge(parent, child, getProbabilityText(child));
                    LinkedList<Edge> inEdges = new LinkedList<Edge>();
                    inEdges.add(edge);
                    _networkNodeToInEdges.put(child, inEdges);
                    _edges.add(edge);
                    _networkNodes.add(child);
                    _networkNodeToPrimarySyntaxNode.put(child, child);
                    return null;

                }

                public Object forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input) {
                    forHybridNodeQualifierNotEmpty(qualifier);
                    return null;
                }

                public Object forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input)  {
                    forHybridNodeQualifierNotEmpty(qualifier);
                    return null;
                }

                private void forHybridNodeQualifierNotEmpty(HybridNodeQualifierNonEmpty qualifier) {

                    String hybridNodeIndex = qualifier.HybridNodeIndex.Content;

                    LinkedList<Edge> inEdges = null;
                    if(_networkNodeToInEdges.containsKey(hybridNodeIndex))
                    {
                       inEdges = _networkNodeToInEdges.get(hybridNodeIndex);
                    }
                    else
                    {
                        inEdges = new LinkedList<Edge>();
                        _networkNodeToInEdges.put(hybridNodeIndex, inEdges);
                        _networkNodes.add(hybridNodeIndex);
                        _networkNodeToPrimarySyntaxNode.put(hybridNodeIndex, child);
                    }

                    Edge edge = new Edge(parent, child, getProbabilityText(child));
                    inEdges.add(edge);
                    _edges.add(edge);
                }

            }, null);

            gatherInEdges(child, subTree.Descendants);
        }

    }

    public Iterable<Edge> getAllInEdges(Object node) {
        return _networkNodeToInEdges.get(node);
    }

    public Object getTail(Edge edge) {
        return edge.Tail;
    }

    public String getEdgeProbabilityText(Edge edge) {
        return edge.ProbabilityText;
    }

    public String getNodeLabelText(NetworkInfo node) {
        return node.NodeLabel.execute(new NodeLabelAlgo<String, Object, RuntimeException>() {

            public String forNodeLabelNonEmpty(NodeLabelNonEmpty node, Object input)  {

                return node.Label.Content;
            }

            public String forNodeLabelEmpty(NodeLabelEmpty node, Object input)  {
                return null;
            }
        }, null);
    }

    public int getNodeLabelTextLineNumber(NetworkInfo node) {

        return node.NodeLabel.execute(new NodeLabelAlgo<Integer, Object, RuntimeException>() {
            public Integer forNodeLabelNonEmpty(NodeLabelNonEmpty node, Object input)  {

                return node.Label.LineNumberStart;

            }

            public Integer forNodeLabelEmpty(NodeLabelEmpty node, Object input) {
                return -1;
            }
        }, null).intValue();

    }

    public int getNodeLabelTextColumnNumber(NetworkInfo node) {

        return node.NodeLabel.execute(new NodeLabelAlgo<Integer, Object, RuntimeException>() {
            public Integer forNodeLabelNonEmpty(NodeLabelNonEmpty node, Object input) {

                return node.Label.ColumnNumberStart;

            }

            public Integer forNodeLabelEmpty(NodeLabelEmpty node, Object input) {
                return -1;
            }
        }, null).intValue();
    }

    public String getHybridNodeIndexText(NetworkInfo node) {

        return node.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<String, Object, RuntimeException>() {
            public String forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input) throws RuntimeException {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public String forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input)  {
                return qualifier.HybridNodeIndex.Content;
            }

            public String forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input)  {
                return qualifier.HybridNodeIndex.Content;
            }
        }, null);
    }

    public int getHybridNodeIndexLineNumber(NetworkInfo node) {

        return node.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<Integer, Object, RuntimeException>() {
            public Integer forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input)  {
                return -1;
            }

            public Integer forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input) {
                return qualifier.HybridNodeIndex.LineNumberStart;
            }

            public Integer forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input) {
                return qualifier.HybridNodeIndex.LineNumberStart;
            }
        }, null).intValue();

    }

    public int getHybridNodeIndexColumnNumber(NetworkInfo node) {
        return node.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<Integer, Object, RuntimeException>() {
            public Integer forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input)  {
                return -1;
            }

            public Integer forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input) {
                return qualifier.HybridNodeIndex.ColumnNumberStart;
            }

            public Integer forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input) {
                return qualifier.HybridNodeIndex.ColumnNumberStart;
            }
        }, null).intValue();
    }

    public String getHybridNodeType(NetworkInfo node) {
        return node.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<String, Object, RuntimeException>() {
            public String forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input) {
                return null;
            }

            public String forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input) {
                return null;
            }

            public String forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input) {
                return qualifier.HybridNodeType.Content;
            }
        }, null);
    }

    public int getHybridNodeTypeLineNumber(NetworkInfo node) {
         return node.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<Integer, Object, RuntimeException>() {
            public Integer forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input) {
                return -1;
            }

            public Integer forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input) {
                return -1;
            }

            public Integer forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input) {
                return qualifier.HybridNodeType.LineNumberStart;
            }
        }, null).byteValue();
    }

    public int getHybridNodeTypeColumnNumber(NetworkInfo node) {
         return node.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<Integer, Object, RuntimeException>() {
            public Integer forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input) {
                return -1;
            }

            public Integer forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input) {
                return -1;
            }

            public Integer forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input) {
                return qualifier.HybridNodeType.ColumnNumberStart;
            }
        }, null).byteValue();
    }

    public String getBranchLengthText(NetworkInfo node) {
        return node.BranchLength.execute(new BranchLengthAlgo<String, Object, RuntimeException>() {
            public String forBranchLengthEmpty(BranchLengthEmpty branchLength, Object input) {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public String forBranchLengthNonEmpty(BranchLengthNonEmpty branchLength, Object input) {
                return branchLength.Length.Content;
            }
        }, null);
    }

    public int getBranchLengthLineNumber(NetworkInfo node) {
         return node.BranchLength.execute(new BranchLengthAlgo<Integer, Object, RuntimeException>() {
            public Integer forBranchLengthEmpty(BranchLengthEmpty branchLength, Object input) {
                return null;
            }

            public Integer forBranchLengthNonEmpty(BranchLengthNonEmpty branchLength, Object input) {
                return branchLength.Length.LineNumberStart;
            }
        }, null).intValue();
    }

    public int getBranchLengthColumnNumber(NetworkInfo node) {
        return node.BranchLength.execute(new BranchLengthAlgo<Integer, Object, RuntimeException>() {
            public Integer forBranchLengthEmpty(BranchLengthEmpty branchLength, Object input) {
                return null;
            }

            public Integer forBranchLengthNonEmpty(BranchLengthNonEmpty branchLength, Object input) {
                return branchLength.Length.ColumnNumberStart;
            }
        }, null).intValue();
    }

    public String getSupportText(NetworkInfo node) {
        return node.Support.execute(new SupportAlgo<String, Object, RuntimeException>() {
            public String forSupportNonEmpty(SupportNonEmpty support, Object input) {
                return support.SupportValue.Content;
            }

            public String forSupportEmpty(SupportEmpty support, Object input) {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null);
    }

    public int getSupportLineNumber(NetworkInfo node) {
        return node.Support.execute(new SupportAlgo<Integer, Object, RuntimeException>() {
            public Integer forSupportNonEmpty(SupportNonEmpty support, Object input) {
                return support.SupportValue.LineNumberStart;
            }

            public Integer forSupportEmpty(SupportEmpty support, Object input) {
                return -1;  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null).intValue();
    }

    public int getSupportColumnNumber(NetworkInfo node) {
        return node.Support.execute(new SupportAlgo<Integer, Object, RuntimeException>() {
            public Integer forSupportNonEmpty(SupportNonEmpty support, Object input) {
                return support.SupportValue.ColumnNumberStart;
            }

            public Integer forSupportEmpty(SupportEmpty support, Object input) {
                return -1;  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null).intValue();
    }

    public String getProbabilityText(NetworkInfo node) {
        return node.Probability.execute(new ProbabilityAlgo<String, Object, RuntimeException>() {
            public String forProbabilityEmpty(ProbabilityEmpty prob, Object input) {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public String forProbabilityNonEmpty(ProbabilityNonEmpty prob, Object input) {
                return prob.ProbabilityValue.Content;
            }
        }, null);
    }

    public int getProbabilityLineNumber(NetworkInfo node) {
       return node.Probability.execute(new ProbabilityAlgo<Integer, Object, RuntimeException>() {
            public Integer forProbabilityEmpty(ProbabilityEmpty prob, Object input) {
                return -1;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public Integer forProbabilityNonEmpty(ProbabilityNonEmpty prob, Object input) {
                return prob.ProbabilityValue.LineNumberStart;
            }
        }, null);
    }

    public int getProbabilityColumnNumber(NetworkInfo node) {
        return node.Probability.execute(new ProbabilityAlgo<Integer, Object, RuntimeException>() {
            public Integer forProbabilityEmpty(ProbabilityEmpty prob, Object input) {
                return -1;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public Integer forProbabilityNonEmpty(ProbabilityNonEmpty prob, Object input) {
                return prob.ProbabilityValue.ColumnNumberStart;
            }
        }, null);
    }

    public NetworkInfo getPrimarySyntaxNode(Object networkNode)
    {
        return _networkNodeToPrimarySyntaxNode.get(networkNode);
    }


}
