package edu.rice.cs.bioinfo.library.phylogenetics;


import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func5;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import junit.framework.Assert;
import org.junit.Test;

import java.util.HashSet;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/25/12
 * Time: 1:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class MULTreeFactoryTest
{
    @Test
    public void testMakeMULTree()
    {


        PhyloGraph<String> network = new DirectedPhyloGraphDefault<String>();
        network.addNode("A");
        network.addNode("B");
        network.addNode("C");
        network.addNode("D");
        network.addNode("E");
        network.addNode("H");
        network.addNode("Z");
        network.addEdge(new PhyloEdge<String>("A", "B"));
        network.addEdge(new PhyloEdge<String>("A", "C"));
        network.addEdge(new PhyloEdge<String>("B", "D"));
        network.addEdge(new PhyloEdge<String>("C", "E"));
        network.addEdge(new PhyloEdge<String>("B", "H"));
        network.addEdge(new PhyloEdge<String>("C", "H"));
        network.addEdge(new PhyloEdge<String>("H", "Z"));

        Func<Graph<MULTreeFactory.MULTreeNode<String>,PhyloEdge<MULTreeFactory.MULTreeNode<String>>>> makeEmptyGraph = new Func<Graph<MULTreeFactory.MULTreeNode<String>, PhyloEdge<MULTreeFactory.MULTreeNode<String>>>>() {
            public Graph<MULTreeFactory.MULTreeNode<String>, PhyloEdge<MULTreeFactory.MULTreeNode<String>>> execute() {

                return new DirectedPhyloGraphDefault<MULTreeFactory.MULTreeNode<String>>();
            }
        };
        Func5<GraphReadOnly<String,PhyloEdge<String>>, PhyloEdge<String>, Graph<MULTreeFactory.MULTreeNode<String>, PhyloEdge<MULTreeFactory.MULTreeNode<String>>>,
                            MULTreeFactory.MULTreeNode<String>, MULTreeFactory.MULTreeNode<String>, PhyloEdge<MULTreeFactory.MULTreeNode<String>>> makeEdge =
                new Func5<GraphReadOnly<String, PhyloEdge<String>>, PhyloEdge<String>, Graph<MULTreeFactory.MULTreeNode<String>, PhyloEdge<MULTreeFactory.MULTreeNode<String>>>, MULTreeFactory.MULTreeNode<String>, MULTreeFactory.MULTreeNode<String>, PhyloEdge<MULTreeFactory.MULTreeNode<String>>>() {
                    public PhyloEdge<MULTreeFactory.MULTreeNode<String>> execute(GraphReadOnly<String, PhyloEdge<String>> arg1, PhyloEdge<String> arg2, Graph<MULTreeFactory.MULTreeNode<String>, PhyloEdge<MULTreeFactory.MULTreeNode<String>>> arg3, MULTreeFactory.MULTreeNode<String> source, MULTreeFactory.MULTreeNode<String> dest) {
                        return new PhyloEdge<MULTreeFactory.MULTreeNode<String>>(source, dest);
                    }
                };

         Graph<MULTreeFactory.MULTreeNode<String>,PhyloEdge<MULTreeFactory.MULTreeNode<String>>> mulTree =
                 new MULTreeFactory<String, PhyloEdge<String>,PhyloEdge<MULTreeFactory.MULTreeNode<String>>>().makeMULTree(network, makeEmptyGraph, makeEdge);

        Assert.assertEquals(9, IterableHelp.countInt(mulTree.getNodes()));
        Assert.assertEquals(8, IterableHelp.countInt(mulTree.getEdges()));
        HashSet<MULTreeFactory.MULTreeNode<String>> foundNodes = new HashSet<MULTreeFactory.MULTreeNode<String>>(
                IterableHelp.<MULTreeFactory.MULTreeNode<String>,MULTreeFactory.MULTreeNode<String>>toList(mulTree.getNodes()));

        Assert.assertEquals(1, countNodeOccurences(mulTree.getNodes(), "A"));
        Assert.assertEquals(1, countNodeOccurences(mulTree.getNodes(), "B"));
        Assert.assertEquals(1, countNodeOccurences(mulTree.getNodes(), "C"));
        Assert.assertEquals(1, countNodeOccurences(mulTree.getNodes(), "D"));
        Assert.assertEquals(1, countNodeOccurences(mulTree.getNodes(), "E"));
        Assert.assertEquals(2, countNodeOccurences(mulTree.getNodes(), "H"));
        Assert.assertEquals(2, countNodeOccurences(mulTree.getNodes(), "Z"));
        Assert.assertEquals(1, countEdgeFrequency(mulTree, "A", "B"));
        Assert.assertEquals(1, countEdgeFrequency(mulTree, "A", "C"));
        Assert.assertEquals(1, countEdgeFrequency(mulTree, "C", "E"));
        Assert.assertEquals(1, countEdgeFrequency(mulTree, "B", "D"));
        Assert.assertEquals(1, countEdgeFrequency(mulTree, "C", "H"));
        Assert.assertEquals(1, countEdgeFrequency(mulTree, "B", "H"));
        Assert.assertEquals(2, countEdgeFrequency(mulTree, "H", "Z"));

    }

    private int countEdgeFrequency(Graph<MULTreeFactory.MULTreeNode<String>,PhyloEdge<MULTreeFactory.MULTreeNode<String>>> mulTree, String source, String destination)
    {
        int edgeFrequency = 0;
        for(PhyloEdge<MULTreeFactory.MULTreeNode<String>> edge : mulTree.getEdges())
        {
            if(edge.Source.Content.equals(source) && edge.Destination.Content.equals(destination))
                edgeFrequency++;
        }
        return edgeFrequency;
    }

    private int countNodeOccurences(Iterable<MULTreeFactory.MULTreeNode<String>> treeNodes, String value)
    {
        int count = 0;
        for(MULTreeFactory.MULTreeNode<String> treeNode : treeNodes)
        {
            if(treeNode.Content.equals(value))
                count++;
        }

        return count;
    }
}
