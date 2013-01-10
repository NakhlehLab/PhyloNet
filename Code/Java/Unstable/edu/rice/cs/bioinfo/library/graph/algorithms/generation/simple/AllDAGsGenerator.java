package edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.HashSet;
import java.util.LinkedList;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/9/13
 * Time: 2:53 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class AllDAGsGenerator<N>
{
    Iterable<N> _nodes;

    private class EdgeRecord
    {
        public final Tuple<N,N> Edge;
        public boolean Include = false;


        private EdgeRecord(Tuple<N, N> edge) {
            Edge = edge;
        }
    }

    AllDAGsGenerator(Iterable<N> nodes)
    {
        _nodes = nodes;
    }

    public void generate()
    {
        LinkedList<N> nodes = new LinkedList<N>(IterableHelp.toList(_nodes));

        LinkedList<EdgeRecord> edgePermutation = new LinkedList<EdgeRecord>();
        for(N source : nodes)
        {
            for(N destination : nodes)
            {
                if(source != destination)
                {
                    edgePermutation.add(new EdgeRecord(new Tuple<N,N>(source, destination)));
                }
            }
        }

        boolean moreToGenerate = true;
        while(moreToGenerate)
        {
            int numTruesChangedToFalse = 0;
            for(EdgeRecord record : edgePermutation)
            {
                if(!record.Include)
                {
                    record.Include = true;
                    break;
                }
                else
                {
                    record.Include = false;
                    numTruesChangedToFalse++;
                }
            }
            generateGraph(edgePermutation);

            moreToGenerate = numTruesChangedToFalse != edgePermutation.size();
        }


    }

    private void generateGraph(LinkedList<EdgeRecord> edgePermutation)
    {
        LinkedList<Tuple<N,N>> edgesOfGraph = new LinkedList<Tuple<N, N>>();

        for(EdgeRecord r : edgePermutation)
        {
            if(r.Include)
            {
                edgesOfGraph.add(r.Edge);
            }
        }

        graphGenerated(edgesOfGraph);
    }

    protected abstract void graphGenerated(Iterable<Tuple<N,N>> edgesOfGraph);
}
