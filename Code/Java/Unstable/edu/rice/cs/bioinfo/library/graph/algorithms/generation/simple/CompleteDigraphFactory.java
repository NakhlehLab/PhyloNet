package edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/11/13
 * Time: 7:21 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class CompleteDigraphFactory<N,E>
{
    public Set<E> makeCompleteDigraph(Set<N> nodes)
    {
        final HashSet<E> graphEdges = new HashSet<E>();
        for(N source : nodes)
        {
            for(N destination : nodes)
            {
                if(source != destination)
                {
                    graphEdges.add(makeEdge(source,destination));
                }
            }
        }

        return graphEdges;
    }

    public abstract E makeEdge(N source, N destination);
}
