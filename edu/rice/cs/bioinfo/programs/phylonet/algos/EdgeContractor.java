package edu.rice.cs.bioinfo.programs.phylonet.algos;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.LinkedList;

public class EdgeContractor
{
    public static <T, N, E> boolean smoothInlineNodeIfNeeded(T tree, boolean isTreeRooted, Func3<T, N, E, Boolean> isDestinationNode, N node, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N, N>> getNodesOfEdge, Proc2<T, E> removeEdge, Proc2<T, E> addEdge, Func3<T, N, N, E> makeEdge, /*out*/ Ref<E> addedEdge)
    {
        N forget1 = null;
        E forget2 = null;
        Ref<N> referenceToForget1 = new Ref<N>(forget1);
        Ref<E> referenceToForget2 = new Ref<E>(forget2);
        return smoothInlineNodeIfNeeded(tree, isTreeRooted, isDestinationNode, node, getIncidentEdges, getNodesOfEdge, removeEdge, addEdge, makeEdge, /*out*/ referenceToForget1, /*out*/ referenceToForget1, /*out*/ addedEdge, /*out*/ referenceToForget2, /*out*/ referenceToForget2);
    }

    public static <T, N, E> boolean smoothInlineNodeIfNeeded(T tree, boolean isTreeRooted, Func3<T, N, E, Boolean> isDestinationNode, N node, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N, N>> getNodesOfEdge, Proc2<T, E> removeEdge,
                                                           Proc2<T, E> addEdge, Func3<T,N,N,E> makeEdge, /*out*/ Ref<N> smoothedNode1, /*out*/ Ref<N> smoothedNode2, /*out*/ Ref<E> addedEdge, /*out*/ Ref<E> removedEdge1, /*out*/ Ref<E> removedEdge2)
    {
        Object[] adjEdges = IterableHelp.toArray(getIncidentEdges.execute(tree, node));


        if (adjEdges.length == 2)
        {
            E removedEdge1ToBe = (E)adjEdges[0];
            E removedEdge2ToBe = (E)adjEdges[1];

            Tuple<N, N> firstEdgeNodes = getNodesOfEdge.execute(tree, removedEdge1ToBe);
            Tuple<N, N> secondEdgeNodes = getNodesOfEdge.execute(tree, removedEdge2ToBe);
            N smoothedNode1ToBe = firstEdgeNodes.Item1.equals(node)  ? firstEdgeNodes.Item2 : firstEdgeNodes.Item1;
            N smoothedNode2ToBe = secondEdgeNodes.Item1.equals(node) ? secondEdgeNodes.Item2 : secondEdgeNodes.Item1;

            removeEdge.execute(tree, removedEdge1ToBe);
            removeEdge.execute(tree, removedEdge2ToBe);

            if(isTreeRooted && isDestinationNode.execute(tree, smoothedNode1ToBe, removedEdge1ToBe))
            {
               addedEdge.set(makeEdge.execute(tree, smoothedNode2ToBe, smoothedNode1ToBe));
            }
            else
            {
               addedEdge.set(makeEdge.execute(tree, smoothedNode1ToBe, smoothedNode2ToBe));
            }

            addEdge.execute(tree, addedEdge.get());
            removedEdge1.set(removedEdge1ToBe);
            removedEdge2.set(removedEdge2ToBe);
            smoothedNode1.set(smoothedNode1ToBe);
            smoothedNode2.set(smoothedNode2ToBe);
            return true;

        }

        return false;
    }

    public static <T, N, E> void undoSmoothing(T tree, N removedNode, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N, N>> getNodesOfEdge, Proc2<T, E> removeEdge,
                                                Proc2<T, E> addEdge, Func3<T, N, N, E> makeEdge, N smoothedNode1, N smoothedNode2, E removedEdge1, E removedEdge2)
    {
        LinkedList<E> incidentEdges = new LinkedList<E>();

        for(E edge : getIncidentEdges.execute(tree, smoothedNode1))
        {
            incidentEdges.add(edge);
        }

        for (E edge :incidentEdges) // linked list to avoid potential concurrent modification exception
        {
            Tuple<N, N> edgeNodes = getNodesOfEdge.execute(tree, edge);

            if (edgeNodes.Item1.equals(smoothedNode2) || edgeNodes.Item2.equals(smoothedNode2))
            {
                removeEdge.execute(tree, edge);
                addEdge.execute(tree, removedEdge1);
                addEdge.execute(tree, removedEdge2);

            }
        }
    }
}