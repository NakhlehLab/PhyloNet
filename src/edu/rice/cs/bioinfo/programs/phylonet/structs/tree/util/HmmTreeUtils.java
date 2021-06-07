package edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.StringReader;
import java.util.Iterator;

/**
 * Created by ethan_000 on 7/29/2014.
 */
public class HmmTreeUtils
{
    public static Tree toTree(String s)
    {
        NewickReader a = new NewickReader(new StringReader(s));
        try {
            return a.readTree();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public static Tree flattenTree(Tree source)
    {

        STITree copy = new STITree();
        flattenNode(source.getRoot(),copy.getRoot(),0);

        return copy;
    }

    private static void processChildren(Iterator<? extends TNode> children, STINode destination,double lengthCount)
    {
        while(children.hasNext())
        {
            TNode child = children.next();
            if (child.getChildCount() == 1)
                processChildren(child.getChildren().iterator(),destination,lengthCount + child.getParentDistance());
            else if (child.getChildCount() == 2 || child.getChildCount() == 0)
                flattenNode(child,destination.createChild(),lengthCount + child.getParentDistance());
            else
            {
                STINode replacement = destination.createChild();
                replacement.setParentDistance(child.getParentDistance());

                Iterator<? extends TNode> theChildren = child.getChildren().iterator();

                TNode firstChild = theChildren.next();

                flattenNode(firstChild,replacement.createChild(),lengthCount + firstChild.getParentDistance());

                STINode replacementChild = replacement.createChild();
                replacementChild.setParentDistance(0);
                processChildren(theChildren,replacementChild,lengthCount);
            }

        }


    }

    private static void flattenNode(TNode root, STINode destination, double length)
    {

        destination.setName(root.getName());
        destination.setParentDistance(length);

        processChildren(root.getChildren().iterator(),destination,0);

    }

    /**
     * Applies a given array of gene tree lengths to a tree.
     * @param tree The tree whose distances need to be set.
     * @param geneTreeLengths The distances to ues.
     */
    public static void setParentDistances(Tree tree,double[] geneTreeLengths) {
        int i = 0;
        for (TNode n : tree.getNodes())
        {
            if (!n.isRoot())
                n.setParentDistance(geneTreeLengths[i++]);
        }
    }

    /**
     * Determines whether or not two trees have the same topology.
     * @param firstTree The first tree.
     * @param secondTree The second tree.
     * @return True if the trees have the same topology.
     */
    public static boolean sameTopology(Tree firstTree,Tree secondTree)
    {
        return sameTopology(firstTree.getRoot(),secondTree.getRoot());
    }

    /**
     * Checks if the parent has a child equivalent to node.
     * @return True if the parent has an equivalent child to node.
     */
    private static boolean parentHasChild(TNode node, TNode parent)
    {
        for (TNode child : parent.getChildren())
        {
            if (sameTopology(node,child))
                return true;
        }

        return false;
    }

    /**
     * Checks whether two tree nodes have the same topology recursively.
     * @param node1 The first node
     * @param node2 The second node
     * @return True if they nodes have the same.
     */
    private static boolean sameTopology(TNode node1, TNode node2)
    {
        if (node1.isLeaf() && node2.isLeaf())
            return node1.getName().equals(node2.getName());
        else if (!node1.isLeaf() && !node2.isLeaf() && node1.getChildCount() == node2.getChildCount())
        {
            for (TNode child: node1.getChildren())
            {
                if (!parentHasChild(child,node2))
                    return false;
            }
            return true;
        }
        else
            return false;
    }

    public static double getHeightSquaredDifference(Tree tree1, Tree tree2)
    {
        return getHeightSquaredDifference(annotateHeights(tree1),annotateHeights(tree2));
    }

    public static double getHeightSquaredDifference(STITree<Double> tree1, STITree<Double> tree2)
    {
        return getHeightSquaredDifferenceChildren(tree1.getRoot(), tree2.getRoot());
    }

    private static double getHeightSquaredDifferenceChildren(STINode<Double> node1, STINode<Double> node2)
    {
        double sum = 0;
        for (STINode<Double> child1 : node1.getChildren())
        {
            sum += getHeightSquaredDifferenceChild(child1,node2);
        }
        return sum;
    }

    private static double getHeightSquaredDifference(STINode<Double> node1, STINode<Double> node2)
    {
        double currentDiff = (node1.getParentDistance() - node2.getParentDistance()) * (node1.getParentDistance() - node2.getParentDistance());
        return currentDiff + getHeightSquaredDifferenceChildren(node1,node2);
    }

    private static double getHeightSquaredDifferenceChild(STINode<Double> child1,STINode<Double> node2)
    {
        for (STINode<Double> child2: node2.getChildren())
        {
            if (sameTopology(child2,child1))
                return getHeightSquaredDifference(child1,child2);
        }
        throw new RuntimeException("Tree 1 and Tree 2 have different topologies?");
    }

    private static STITree<Double> annotateHeights(Tree tree)
    {
        STITree<Double> copy = new STITree<Double>();
        annotateHeightsNode(copy.getRoot(),tree.getRoot());
        return copy;

    }

    private static void annotateHeightsNode(STINode<Double> destination, TNode source)
    {
        destination.setName(source.getName());
        destination.setParentDistance(source.getParentDistance());

        if (source.isLeaf())
        {
            destination.setData(0.0);
        }
        else
        {
            double myHeight = Double.NaN;
            for (TNode sourceChild : source.getChildren())
            {
                STINode<Double> destinationChild = destination.createChild();
                annotateHeightsNode(destinationChild,sourceChild);
                double myHeightPart = destinationChild.getData() + destinationChild.getParentDistance();

                if (Double.isNaN(myHeight))
                    myHeight = myHeightPart;
                else if (Math.abs(myHeight - myHeightPart) > .1)
                    throw new RuntimeException("Heights should sum correctly");


            }
            destination.setData(myHeight);
        }

    }
}
