package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree.nni;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.IsDestinationNode;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Proc4;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/4/12
 * Time: 4:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class NearestNeighborInterchangeInPlace<T extends Graph<N,E>, N,E> extends NearestNeighborInterchangeBase<T, N,E>
{
    public NearestNeighborInterchangeInPlace(Func2<N, N, E> makeEdge) {
        super(makeEdge);
    }

    @Override
    public void computeRearrangementsWithoutValidation(T tree, Proc4<T,E,E,E> rearrangementComputed) {

        LinkedList<E> edges = new LinkedList<E>(); // get a list to prevent concurrent modification exception
        for(E edge : tree.getEdges())
        {
            edges.add(edge);
        }

        for (E edge : edges)
        {
            Tuple<N,N> nodesOfEdge = tree.getNodesOfEdge(edge);
            int node1Degree = IterableHelp.countInt(tree.getIncidentEdges(nodesOfEdge.Item1));
            int node2Degree = IterableHelp.countInt(tree.getIncidentEdges(nodesOfEdge.Item2));
            if (node1Degree > 1 && node2Degree > 1) // is internal edge
            {
                if (node1Degree == 3 && node2Degree == 3)
                {
                    swapDoubleDegreeThree(tree, rearrangementComputed, edge, nodesOfEdge.Item1, nodesOfEdge.Item2);
                }
                else
                {
                    swapGeneralCase(tree, rearrangementComputed, edge, nodesOfEdge.Item1, nodesOfEdge.Item2);
                }
            }
        }
    }

    private void swapGeneralCase(T tree, Proc4<T,E,E,E> rearrangementComputed, E internalEdge, N internalEdgeNodeA, N internalEdgeNodeB)
    {
        LinkedList<E> swapEdgesA = new LinkedList<E>();

        for(E edge : tree.getIncidentEdges(internalEdgeNodeA))
        {
            if(!edge.equals(internalEdge))
            {
                swapEdgesA.add(edge);
            }
        }

        LinkedList<E> swapEdgesB = new LinkedList<E>();

        for(E edge : tree.getIncidentEdges(internalEdgeNodeB))
        {
            if(!edge.equals(internalEdge))
            {
                swapEdgesB.add(edge);
            }
        }

        swapHelp(tree, rearrangementComputed, internalEdge, swapEdgesA, internalEdgeNodeA, swapEdgesB, internalEdgeNodeB, true, true);
    }

    private void swapDoubleDegreeThree(T tree, Proc4<T,E,E,E> rearrangementComputed, E internalEdge, N internalEdgeNodeA, N internalEdgeNodeB) {

        LinkedList<E> swapEdgesA = new LinkedList<E>();
        IsDestinationNode isDestinationNode = new IsDestinationNode();
        for(E edge : tree.getIncidentEdges(internalEdgeNodeA))
        {
            if(!edge.equals(internalEdge))
            {
                if(tree.isRooted())
                {
                    if(!isDestinationNode.execute(tree, internalEdgeNodeA, edge))
                    {
                        swapEdgesA.add(edge);
                    }
                }
                else
                {
                    swapEdgesA.add(edge);
                }
            }
        }

        LinkedList<E> swapEdgesB = new LinkedList<E>();

        for(E edge : tree.getIncidentEdges(internalEdgeNodeB))
        {
            if(!edge.equals(internalEdge))
            {
                if(tree.isRooted())
                {
                    if(!isDestinationNode.execute(tree, internalEdgeNodeB, edge))
                    {
                       swapEdgesB.add(edge);
                    }
                }
                else
                {
                    swapEdgesB.add(edge);
                }
            }

            if(swapEdgesA.size() == 2 && swapEdgesB.size() == 1)
            {
                break;
            }
        }

        swapHelp(tree, rearrangementComputed, internalEdge, swapEdgesA, internalEdgeNodeA, swapEdgesB, internalEdgeNodeB, true, true);
    }

    private void swapHelp(T tree, Proc4<T,E,E,E> rearrangementComputed, E internalEdge, Iterable<E> swapEdgesA, N internalEdgeNodeA, Iterable<E> swapEdgesB, N internalEdgeNodeB,
                          boolean performUndo, boolean invokeRearrangementComputed)
    {
        for (E swapA : swapEdgesA)
        {
            Tuple<N,N> nodesOfSwapA = tree.getNodesOfEdge(swapA);
            N transplantInA = nodesOfSwapA.Item1.equals(internalEdgeNodeA) ? nodesOfSwapA.Item2 : nodesOfSwapA.Item1;
            tree.removeEdge(swapA);
            E newEdge1 = makeEdge.execute(internalEdgeNodeB, transplantInA);
            tree.addEdge(newEdge1);
            for (E swapB : swapEdgesB)
            {
                Tuple<N, N> nodesOfSwapB = tree.getNodesOfEdge(swapB);
                N transplantInB = nodesOfSwapB.Item1.equals(internalEdgeNodeB) ? nodesOfSwapB.Item2 : nodesOfSwapB.Item1;

                tree.removeEdge(swapB);
                E newEdge2 = makeEdge.execute(internalEdgeNodeA, transplantInB);
                tree.addEdge(newEdge2);

                if(invokeRearrangementComputed)
                    rearrangementComputed.execute(tree, internalEdge, swapA, swapB);

                if(performUndo)
                {
                    tree.removeEdge(newEdge2);
                    tree.addEdge(swapB);
                }
            }

            if(performUndo)
            {
                tree.removeEdge(newEdge1);
                tree.addEdge(swapA);
            }
        }
    }

    public T performInterchange(T tree, E internalEdge, E swapEdgeA, E swapEdgeB)
    {
        Tuple<N,N> nodesOfInternalEdge = tree.getNodesOfEdge(internalEdge);
        Tuple<N,N> nodesOfSwapA = tree.getNodesOfEdge(swapEdgeA);
        Tuple<N,N> nodesOfSwapB = tree.getNodesOfEdge(swapEdgeB);

        N internalEdgeNodeA = findCommonNode(tree, internalEdge, swapEdgeA);
        N internalEdgeNodeB = findCommonNode(tree, internalEdge, swapEdgeB);

        List<E> swapEdgesA = new ArrayList<E>();
        List<E> swapEdgesB = new ArrayList<E>();

        swapEdgesA.add(swapEdgeA);
        swapEdgesB.add(swapEdgeB);

        swapHelp(tree, null, internalEdge, swapEdgesA, internalEdgeNodeA, swapEdgesB, internalEdgeNodeB, false, false);


       return tree;
    }

    private N findCommonNode(T tree, E internalEdge, E swapEdge)
    {
        Tuple<N,N> nodesOfInternalEdge = tree.getNodesOfEdge(internalEdge);
        Tuple<N,N> nodesOfSwap = tree.getNodesOfEdge(swapEdge);

        if(nodesOfSwap.Item1 == nodesOfInternalEdge.Item1)
        {
            return nodesOfInternalEdge.Item1;
        }
        else if(nodesOfSwap.Item2 == nodesOfInternalEdge.Item1)
        {
            return nodesOfInternalEdge.Item1;
        }
         if(nodesOfSwap.Item1 == nodesOfInternalEdge.Item2)
        {
            return nodesOfInternalEdge.Item2;
        }
        else if(nodesOfSwap.Item2 == nodesOfInternalEdge.Item2)
        {
            return nodesOfInternalEdge.Item2;
        }
        else
        {
            throw new IllegalArgumentException("Internal edge and swapEdge do not share an incident node.");
        }
    }

}