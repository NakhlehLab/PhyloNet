package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.*;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.NetNodes;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;

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

            return new NetworkNonEmpty(new RootageQualifierNonEmpty(RootageQualifier.ROOTED),
                                       new DescendantList(subtrees),
                                       rootInfo);
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
        return network.execute(new NetworkAlgo<STITree<T>, Object, RuntimeException>() {

            public STITree<T> forNetworkEmpty(NetworkEmpty networkEmpty, Object o) throws RuntimeException {

                return new STITree<T>();

            }

            public STITree<T> forNetworkNonEmpty(NetworkNonEmpty network, Object o) throws RuntimeException {

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
        }, null);







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
}
