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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkAlgo;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/9/11
 * Time: 12:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkTransformer {


    static <T> Network fromClassicNetwork(edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network<T> classicNetwork)
    {
        if(classicNetwork.isEmpty())
        {
            return NetworkEmpty.Singleton;
        }
        else
        {
            NetNode<T> root = classicNetwork.getRoot();
            ArrayList<NetNode<T>> nodeToHybridIndex = new ArrayList<NetNode<T>>();
            NetworkInfo rootInfo = makeInfo(root, null, nodeToHybridIndex);
            Iterable<Subtree> subtrees = makeSubtrees(root, nodeToHybridIndex);

            return new NetworkNonEmpty(new RootageQualifierNonEmpty(RootageQualifier.ROOTED), new DescendantList(subtrees), rootInfo, TreeProbabilityEmpty.Singleton
            );
        }
    }

    private static <T> Iterable<Subtree> makeSubtrees(NetNode<T> root, ArrayList<NetNode<T>> nodeToHybridIndex) {

        LinkedList<Subtree> tbr = new LinkedList<Subtree>();

        for(NetNode<T> child : root.getChildren())
        {
            NetworkInfo childInfo = makeInfo(child, root, nodeToHybridIndex);
            DescendantList dl = new DescendantList(makeSubtrees(child, nodeToHybridIndex));

            tbr.add(new Subtree(dl, childInfo));
        }

        return tbr;
    }

    private static <T> NetworkInfo makeInfo(NetNode<T> node, NetNode<T> parent, ArrayList<NetNode<T>> nodeToHybridIndex)
    {
        boolean isHybrid = node.getIndeg() > 1;
        int hybridIndex = -1;

        HybridNodeQualifier hybridQuualifier = null;
        if(isHybrid)
        {
            if(nodeToHybridIndex.contains(node))
            {
                hybridIndex = nodeToHybridIndex.indexOf(node) + 1;
            }
            else
            {
                hybridIndex = nodeToHybridIndex.size() + 1;
                nodeToHybridIndex.add(node);
            }

            hybridQuualifier = new HybridNodeQualifierNonEmpty(new Text(hybridIndex + "", -1, -1, false));
        }
        else
        {
            hybridQuualifier = HybridNodeQualifierEmpty.Singleton;
        }

        BranchLength bl = BranchLengthEmpty.Singleton;

        if(parent != null)
        {
            double parentDistance = node.getParentDistance(parent);

            if(parentDistance > 0)
            {
                bl = new BranchLengthNonEmpty(new Text(parentDistance + "", -1, -1, false));
            }

        }

        return new NetworkInfo(new NodeLabelNonEmpty(new Text(node.getName(), -1, -1, false)),
                               hybridQuualifier, bl , SupportEmpty.Singleton, ProbabilityEmpty.Singleton );
    }

    static STITree<Object> toTree(Network network)
    {
        return toTree(network, new Func1<NetworkInfo, Object>() {
            public Object execute(NetworkInfo networkInfo) {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }
        });
    }

    static <T> STITree<T> toTree(Network network, final Func1<NetworkInfo, T> getData)
    {
        return network.execute(new NetworkAlgo<STITree<T>, RuntimeException>() {

            public STITree<T> forNetworkEmpty(NetworkEmpty networkEmpty) throws RuntimeException {

                return new STITree<T>();

            }

            public STITree<T> forNetworkNonEmpty(NetworkNonEmpty network) throws RuntimeException {

                if(network.PrincipleInfo == null)
                {
                    STITree<T> empytTree = new STITree<T>();
                    return empytTree;
                }

                final STITree<T> tbr =  network.RootageQualifier.execute(new RootageQualifierAlgo<STITree<T>, Object, RuntimeException>() {
                    public STITree<T> forEmptyQualifier(RootageQualifierEmpty rootageQualifierEmpty, Object o) throws RuntimeException {
                        return new STITree<T>(true);
                    }

                    public STITree<T> forNonEmptyQualifier(RootageQualifierNonEmpty rootageQualifierNonEmpty, Object o) throws RuntimeException {

                        return new STITree<T>(rootageQualifierNonEmpty.isRooted());
                    }
                }, null);

                STINode<T> primaryNode = tbr.getRoot();
                setNodeInfo(primaryNode, network.PrincipleInfo);
                primaryNode.setData(getData.execute(network.PrincipleInfo));

                HashSet<NetworkInfo> seenNodes = new HashSet<NetworkInfo>();
                seenNodes.add(network.PrincipleInfo);

                buildTree(primaryNode, network.PrincipleDescendants, seenNodes, getData);

                return tbr;
            }
        });

    }

    private static <T> void buildTree(STINode<T> parent, DescendantList children,
                                  HashSet<NetworkInfo> seenNodes, Func1<NetworkInfo, T> getData)
    {

        for(Subtree child : children.Subtrees)
        {
            try
            {
                seenNodes.add(child.NetworkInfo);
            }
            catch(IllegalArgumentException e)
            {
                throw new IllegalArgumentException("Passed network is cyclic.");
            }

            STINode<T> treeChild = parent.createChild();
            setNodeInfo(treeChild, child.NetworkInfo);
            treeChild.setData(getData.execute(child.NetworkInfo));

            buildTree(treeChild, child.Descendants, seenNodes, getData);
        }
    }

    private static <T> void setNodeInfo(final STINode<T> treeNode, NetworkInfo networkNode)
    {
        String nodeName = networkNode.NodeLabel.execute(new NodeLabelAlgo<String, Object, RuntimeException>() {
            public String forNodeLabelNonEmpty(NodeLabelNonEmpty nodeLabelNonEmpty, Object o) throws RuntimeException {
                return nodeLabelNonEmpty.Label.Content;
            }

            public String forNodeLabelEmpty(NodeLabelEmpty nodeLabelEmpty, Object o) throws RuntimeException {
                return STITree.NO_NAME;  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null);

        networkNode.BranchLength.execute(new BranchLengthAlgo<Object, Object, RuntimeException>() {
            public Object forBranchLengthEmpty(BranchLengthEmpty branchLengthEmpty, Object o) throws RuntimeException {

                treeNode.setParentDistance(STINode.NO_DISTANCE);
                return null;

            }

            public Object forBranchLengthNonEmpty(BranchLengthNonEmpty branchLengthNonEmpty, Object o) throws RuntimeException {

                treeNode.setParentDistance(Double.parseDouble(branchLengthNonEmpty.Length.Content));
                return null;

            }
        }, null);

        treeNode.setName(nodeName);
    }

    public static String toENewickTree(NetworkNonEmpty network)
    {
        if(network.execute(ContainsHybridNode.Singleton, null))
        {
            throw new IllegalArgumentException("Passed network is not a tree.");
        }

        SingleLinePrinter printer = new SingleLinePrinter();
        printer.setExcludeProbability(true);
        //printer.setSupportTransformer(new TransformSupportToBase100());
        return printer.toString(network);

    }

    public static String toENewick(NetworkNonEmpty network)
    {
        if(!network.execute(ContainsHybridNode.Singleton, null))
        {
            return toENewickTree(network);
        }

        Map<String, NetworkNonEmpty> hybridIndexToSubNetwork = new HashMap<String, NetworkNonEmpty>();
        buildSubNetworks(hybridIndexToSubNetwork, network.PrincipleInfo, network.PrincipleDescendants);
        Map<String, String> hybridIndexToHybridNodeLabel = makeHybridIndexToHybridNodeLabelMapping(network);

        NetworkNonEmpty relabeledNetwork = makeRelabeledNetwork(network, hybridIndexToHybridNodeLabel);

        SingleLinePrinter printer = new SingleLinePrinter();
        //printer.setSupportTransformer(new TransformSupportToBase100());

        StringBuilder b = new StringBuilder();
        b.append("N = " + printer.toString(relabeledNetwork));

        for(String hybridNodeIndex : hybridIndexToSubNetwork.keySet())
        {
            NetworkNonEmpty subNetwork =  hybridIndexToSubNetwork.get(hybridNodeIndex);
            String rootLabel =  hybridIndexToHybridNodeLabel.get(hybridNodeIndex);
            NetworkNonEmpty relabeledSubNetwork = makeRelabeledNetwork(subNetwork, hybridIndexToHybridNodeLabel);
            String relabeledString = printer.toString(relabeledSubNetwork);

      //      NetworkNonEmpty subNetworkDescList = new NetworkNonEmpty(new RootageQualifierNonEmpty(RootageQualifier.ROOTED) , subNetwork.PrincipleDescendants,
        //            new NetworkInfo(NodeLabelEmpty.Singleton, HybridNodeQualifierEmpty.Singleton, BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton));

            // at this point relabeledString will be a newick style string that looks something like:
            //        (A,B)R;
            // but in eNewick for the repeated trees of a network we need the form:
            //        R = (A,B);
            String relabeledStringWithoutApendedRoot = relabeledString.substring(0, relabeledString.lastIndexOf(')') + 1) + ";";
            b.append("\n" + rootLabel + " = " + relabeledStringWithoutApendedRoot);
        }

        return b.toString();
    }

    private static NetworkNonEmpty makeRelabeledNetwork(NetworkNonEmpty network, Map<String, String> hybridIndexToHybridNodeLabel) {

        return new NetworkNonEmpty(network.RootageQualifier,
                makeRelabeledDL(network.PrincipleDescendants, hybridIndexToHybridNodeLabel),
                makeRelabeledNode(network.PrincipleInfo, hybridIndexToHybridNodeLabel));

    }

    private static DescendantList makeRelabeledDL(DescendantList descendents, Map<String, String> hybridIndexToHybridNodeLabel)
    {
        LinkedList<Subtree> trees = new LinkedList<Subtree>();

        for(Subtree t : descendents.Subtrees)
        {
            Boolean isHybridNode = t.NetworkInfo.HybridNodeQualifier.execute(new IsHybridNode(), null);

            NetworkInfo newInfo = makeRelabeledNode(t.NetworkInfo, hybridIndexToHybridNodeLabel);
            DescendantList newDl = isHybridNode ? DescendantList.EMPTY_DESCENDANT_LIST : makeRelabeledDL(t.Descendants, hybridIndexToHybridNodeLabel);
            trees.add(new Subtree(newDl,newInfo));
        }

        return new DescendantList(trees);
    }

    private static NetworkInfo makeRelabeledNode(final NetworkInfo node, final Map<String, String> hybridIndexToHybridNodeLabel) {

       return node.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<NetworkInfo, Object, RuntimeException>() {

           public NetworkInfo forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty hybridNodeQualifierEmpty, Object o) throws RuntimeException {
               return node;
           }

           public NetworkInfo forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty hybridNodeQualifierNonEmpty, Object o) throws RuntimeException {

               NodeLabel label = new NodeLabelNonEmpty(new Text(hybridIndexToHybridNodeLabel.get(hybridNodeQualifierNonEmpty.HybridNodeIndex.Content), -1, -1, false));
               return new NetworkInfo(label, HybridNodeQualifierEmpty.Singleton, node.BranchLength, node.Support, node.Probability);

           }

           public NetworkInfo forHybridNodeQualifierWithType(HybridNodeQualifierWithType hybridNodeQualifierWithType, Object o) throws RuntimeException {

               return forHybridNodeQualifierNonEmpty(hybridNodeQualifierWithType, o);
           }

       }, null);
    }

    private static Map<String, String> makeHybridIndexToHybridNodeLabelMapping(NetworkNonEmpty network)
    {
        HashMap<String,String> hybridIndexToHybridNodeLabel = new HashMap<String, String>();

         HashSet<String> labelSet = new HashSet<String>();

        for(String label : network.execute(new ExtractNodeLabels(), null))
        {
            if(!labelSet.contains(label))
            {
                labelSet.add(label);
            }
        }

       makeHybridIndexToHybridNodeLabelMappingHelp(network.PrincipleInfo, network.PrincipleDescendants, labelSet, hybridIndexToHybridNodeLabel);

        return hybridIndexToHybridNodeLabel;


    }

    private static void makeHybridIndexToHybridNodeLabelMappingHelp(NetworkInfo node, DescendantList children, final HashSet<String> labelSet,
                                                                    final HashMap<String, String> hybridIndexToHybridNodeLabel)
    {
       final String hybridNodeIndex = node.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<String, Object, RuntimeException>() {
            public String forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty hybridNodeQualifierEmpty, Object o) throws RuntimeException {
                return null;
            }

            public String forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty hybridNodeQualifierNonEmpty, Object o) throws RuntimeException {
                return hybridNodeQualifierNonEmpty.HybridNodeIndex.Content;
            }

            public String forHybridNodeQualifierWithType(HybridNodeQualifierWithType hybridNodeQualifierWithType, Object o) throws RuntimeException {
                 return hybridNodeQualifierWithType.HybridNodeIndex.Content;
            }
        }, null);

        if(hybridNodeIndex != null)
        {
            node.NodeLabel.execute(new NodeLabelAlgo<Object, Object, RuntimeException>() {
                public Object forNodeLabelNonEmpty(NodeLabelNonEmpty nodeLabelNonEmpty, Object o) throws RuntimeException {

                    hybridIndexToHybridNodeLabel.put(hybridNodeIndex, nodeLabelNonEmpty.Label.Content);
                    return null;
                }

                public Object forNodeLabelEmpty(NodeLabelEmpty nodeLabelEmpty, Object o) throws RuntimeException {

                    String label;
                    do
                    {
                        label = UUID.randomUUID().toString().replace("-", "");
                    }
                    while(labelSet.contains(label));
                    hybridIndexToHybridNodeLabel.put(hybridNodeIndex, label);

                    return null;


                }
            }, null);
        }

        for(Subtree t: children.Subtrees)
        {
            makeHybridIndexToHybridNodeLabelMappingHelp(t.NetworkInfo, t.Descendants, labelSet, hybridIndexToHybridNodeLabel);
        }


    }

    private static void buildSubNetworks(final Map<String, NetworkNonEmpty> hybridIndexToSubNetwork, final NetworkInfo node, final DescendantList children)
    {

        node.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<Object, Object, RuntimeException>()
        {
            public Object forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty hybridNodeQualifierEmpty, Object o) throws RuntimeException {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public Object forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty hybridNodeQualifierNonEmpty, Object o) throws RuntimeException {

                String hybridNodeIndex = hybridNodeQualifierNonEmpty.HybridNodeIndex.Content;

                if(!hybridIndexToSubNetwork.containsKey(hybridNodeIndex))
                {
                    NetworkNonEmpty subNetwork = new NetworkNonEmpty(new RootageQualifierNonEmpty(RootageQualifier.ROOTED), children, node);
                    hybridIndexToSubNetwork.put(hybridNodeIndex, subNetwork);
                }
                else if(children.Subtrees.iterator().hasNext())
                {
                     NetworkNonEmpty subNetwork = new NetworkNonEmpty(new RootageQualifierNonEmpty(RootageQualifier.ROOTED), children, node);
                     hybridIndexToSubNetwork.put(hybridNodeIndex, subNetwork);
                }

                return null;
            }

            public Object forHybridNodeQualifierWithType(HybridNodeQualifierWithType hybridNodeQualifierWithType, Object o) throws RuntimeException {
                return forHybridNodeQualifierNonEmpty(hybridNodeQualifierWithType, o);
            }
        }, null);

        for(Subtree tree : children.Subtrees)
        {
            buildSubNetworks(hybridIndexToSubNetwork, tree.NetworkInfo, tree.Descendants);
        }
    }




}
