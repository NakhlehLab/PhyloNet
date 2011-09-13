package edu.rice.bioinfo.programs.phylonet.commands;

import edu.rice.bioinfo.library.language.richnewick._1_0.ast.*;
import edu.rice.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.HashSet;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/9/11
 * Time: 12:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkTransformer {

    static Tree toTree(Network network)
    {
        if(network.PrincipleInfo == null)
        {
            STITree<Object> empytTree = new STITree<Object>();
            return empytTree;
        }

        final STITree<Object> tbr =  network.RootageQualifier.execute(new RootageQualifierAlgo<STITree<Object>, Object, RuntimeException>() {
            public STITree<Object> forEmptyQualifier(RootageQualifierEmpty rootageQualifierEmpty, Object o) throws RuntimeException {
               return new STITree<Object>(true);
            }

            public STITree<Object> forNonEmptyQualifier(RootageQualifierNonEmpty rootageQualifierNonEmpty, Object o) throws RuntimeException {

                return new STITree<Object>(rootageQualifierNonEmpty.isRooted());
            }
        }, null);

        STINode<Object> primaryNode = tbr.getRoot();
        setNodeName(primaryNode, network.PrincipleInfo);




        HashSet<NetworkInfo> seenNodes = new HashSet<NetworkInfo>();
        seenNodes.add(network.PrincipleInfo);

        buildTree(primaryNode, network.PrincipleDescendants, seenNodes);

        return tbr;
    }

    private static void buildTree(STINode<Object> parent, DescendantList children,
                                  HashSet<NetworkInfo> seenNodes)
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

            STINode<Object> treeChild = parent.createChild();
            setNodeName(treeChild, child.NetworkInfo);

            buildTree(treeChild, child.Descendants, seenNodes);
        }
    }

    private static void setNodeName(STINode<Object> treeNode, NetworkInfo networkNode)
    {
           String nodeName = networkNode.NodeLabel.execute(new NodeLabelAlgo<String, Object, RuntimeException>() {
            public String forNodeLabelNonEmpty(NodeLabelNonEmpty nodeLabelNonEmpty, Object o) throws RuntimeException {
                return nodeLabelNonEmpty.Label.Content;
            }

            public String forNodeLabelEmpty(NodeLabelEmpty nodeLabelEmpty, Object o) throws RuntimeException {
                return STITree.NO_NAME;  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null);

        treeNode.setName(nodeName);
    }
}
