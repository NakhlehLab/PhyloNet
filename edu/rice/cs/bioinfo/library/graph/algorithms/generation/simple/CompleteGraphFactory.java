package edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/12/13
 * Time: 6:02 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class CompleteGraphFactory<N,E>
{

    public Set<E> makeCompleteGraph(Set<N> nodes)
    {
        LinkedList<N> nodeList = new LinkedList<N>(nodes);

        final HashSet<E> graphEdges = new HashSet<E>();
        int node1Index = 0;
        for(N node1 : nodeList)
        {
            ListIterator<N> node2Elements = nodeList.listIterator(node1Index+1);

            while(node2Elements.hasNext())
            {
                N node2 = node2Elements.next();
                graphEdges.add(makeEdge(node1,node2));
            }
            node1Index++;
        }

        return graphEdges;
    }

    public abstract E makeEdge(N node1, N node2);

}
