package edu.rice.cs.bioinfo.library.phylogenetics;


import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/5/12
 * Time: 5:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class AreSameTopology
{
    public class Result
    {

        public final boolean SameTopology;

        public final GraphReadOnly OffendingNetwork;

        public final Object OffendingEdge;

        Result()
        {
            SameTopology = true;
            OffendingNetwork = null;
            OffendingEdge = null;

        }

        Result(GraphReadOnly offendingNetwork, Object offendingEdge) {
            SameTopology = false;
            OffendingNetwork = offendingNetwork;
            OffendingEdge = offendingEdge;
        }
    }

    public <N1,E1,N2> Result execute(GraphReadOnly<N1,E1> graph1, GraphReadOnly graph2, Map<N1,N2> graph1NodeToGraph2Node)
    {
        if(IterableHelp.countInt(graph1.getNodes()) != graph1NodeToGraph2Node.size())
        {
            throw new IllegalArgumentException("Size mismatch between graph1 node set and node map.");
        }

        if(IterableHelp.countInt(graph2.getNodes()) != graph1NodeToGraph2Node.size())
        {
            throw new IllegalArgumentException("Size mismatch between graph2 node set and node map.");
        }

        int graph1EdgeCountAccum = 0;
        Set<Object> matchedEdges = new HashSet<Object>();
        for(E1 graph1Edge : graph1.getEdges())
        {
            graph1EdgeCountAccum++;

            Tuple<N1,N1> nodesOfEdge = graph1.getNodesOfEdge(graph1Edge);

            N2 correspondingSource = graph1NodeToGraph2Node.get(nodesOfEdge.Item1);
            N2 correspondingDest   = graph1NodeToGraph2Node.get(nodesOfEdge.Item2);

            if(graph2.containsEdge(correspondingSource, correspondingDest))
            {
                Object matchedEdge = graph2.getEdge(correspondingSource, correspondingDest);
                matchedEdges.add(matchedEdge);

            }
            else
            {
                return new Result(graph1, graph1Edge);
            }
        }

        int graph1EdgeCount = graph1EdgeCountAccum;

        if(IterableHelp.countInt(graph2.getEdges()) != graph1EdgeCount)
        {
            for(Object graph2Edge : graph2.getEdges())
            {
                if(!matchedEdges.contains(graph2Edge))
                {
                    return new Result(graph2, graph2Edge);
                }
            }
        }

        throw new RuntimeException("Logical Exception:  Should never be reached.");

    }
}
