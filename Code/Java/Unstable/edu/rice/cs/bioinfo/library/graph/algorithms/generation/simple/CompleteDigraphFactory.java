package edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/11/13
 * Time: 7:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class CompleteDigraphFactory<N>
{
    public Set<Tuple<N,N>> makeCompleteDigraph(Set<N> nodes)
    {
        final HashSet<Tuple<N,N>> graphEdges = new HashSet<Tuple<N, N>>();
        for(N source : nodes)
        {
            for(N destination : nodes)
            {
                if(source != destination)
                {
                    graphEdges.add(new Tuple<N,N>(source, destination));
                }
            }
        }

        return graphEdges;
    }
}
