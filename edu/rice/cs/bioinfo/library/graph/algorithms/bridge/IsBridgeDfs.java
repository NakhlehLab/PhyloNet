package edu.rice.cs.bioinfo.library.graph.algorithms.bridge;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.HashSet;
import java.util.Stack;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/10/13
 * Time: 11:32 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class IsBridgeDfs<G,E,N> implements IsBridge<G,E>
{
    public boolean execute(G graph, E edge)
    {
        HashSet<N> seenNodes = new HashSet<N>();
        Stack<N> toExplore = new Stack<N>();

        Tuple<? extends N,? extends N> nodesOfEdge = getNodesOfEdge(edge, graph);
        N searchRoot = nodesOfEdge.Item1;
        toExplore.add(searchRoot);
        seenNodes.add(searchRoot);

        while(!toExplore.isEmpty())
        {
            N node = toExplore.pop();

            for(E incident : getIncidentEdges(node, graph))
            {
                if(edge.equals(incident))
                    continue;

                nodesOfEdge = getNodesOfEdge(incident, graph);
                N neighbor = (N) nodesOfEdge.other(node);

                if(seenNodes.add(neighbor))
                {
                    toExplore.add(neighbor);
                }
            }
        }

        return getNodeCount(graph) != seenNodes.size();
    }

    protected abstract int getNodeCount(G graph);

    protected abstract Iterable<? extends E> getIncidentEdges(N node, G graph);

    protected abstract Tuple<? extends N, ? extends N> getNodesOfEdge(E edge, G graph);
}
