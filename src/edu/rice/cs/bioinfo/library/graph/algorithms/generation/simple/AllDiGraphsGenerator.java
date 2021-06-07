package edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple;

import edu.rice.cs.bioinfo.library.math.discrete.Configurations;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/9/13
 * Time: 2:53 PM
 * To change this template use File | Settings | File Templates.
 */
public class AllDiGraphsGenerator<N> implements Iterable<Set<Tuple<N,N>>>
{
    private Iterable<N> _nodes;

    private class EdgeRecord
    {
        public final Tuple<N,N> Edge;
        public boolean Include = false;


        private EdgeRecord(Tuple<N, N> edge) {
            Edge = edge;
        }
    }

    public AllDiGraphsGenerator(Iterable<N> nodes)
    {
        _nodes = nodes;
    }

    public Iterator<Set<Tuple<N, N>>> iterator()
    {
        final LinkedList<Tuple<N,N>> allPossibleEdges = new LinkedList<Tuple<N, N>>();
        for(N source : _nodes)
        {
            for(N destination : _nodes)
            {
                if(source != destination)
                {
                    allPossibleEdges.add(new Tuple<N,N>(source, destination));
                }
            }
        }


        List<Boolean> inclusionOptions = Arrays.asList(Boolean.TRUE,  Boolean.FALSE);
        LinkedList<List<Boolean>> inclusionFilter = new LinkedList<List<Boolean>>();
        for(Object possibleEdge : allPossibleEdges)
        {
            inclusionFilter.add(inclusionOptions);
        }

       final Iterator<List<Boolean>> inclusionConfigurations = new Configurations<Boolean>(inclusionFilter).iterator();

       return new Iterator<Set<Tuple<N, N>>>() {
           public boolean hasNext() {
               return inclusionConfigurations.hasNext();
           }

           public Set<Tuple<N, N>> next() {
               List<Boolean> edgesToInclude = inclusionConfigurations.next();

               Set<Tuple<N,N>> includedEdges = new HashSet<Tuple<N, N>>();

               Iterator<Tuple<N,N>> edges =  allPossibleEdges.iterator();
               for(Boolean toInclude : edgesToInclude)
               {
                   Tuple<N,N> edge = edges.next();
                   if(toInclude)
                   {
                       includedEdges.add(edge);
                   }

               }

               return includedEdges;
           }

           public void remove() {
               throw new UnsupportedOperationException();
           }
       };

    }


}
