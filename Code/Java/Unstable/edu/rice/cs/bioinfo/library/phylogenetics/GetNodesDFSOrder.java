package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func1;

import java.util.HashSet;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/27/12
 * Time: 12:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class GetNodesDFSOrder<N,E> implements Func1<GraphReadOnly<N,E>, Iterable<N>>
{
    private final GetDirectSuccessors<N,E> _getDirectSuccessors = new GetDirectSuccessors<N, E>();

    public Iterable<N> execute(GraphReadOnly<N,E> graph)
    {
        if(!graph.isRooted())
        {
            throw new IllegalArgumentException("Given graph must be rooted");
        }

        LinkedList<N> nodesDFSOrder = new LinkedList<N>();

        processGraph(graph, new FindRoot<N>().execute(graph), nodesDFSOrder, new HashSet<N>());

        return nodesDFSOrder;
    }

    private void processGraph(GraphReadOnly<N, E> graph, N node, LinkedList<N> nodesDFSOrder, HashSet<N> storedNodes)
    {
        if(!storedNodes.contains(node))
        {
            nodesDFSOrder.addLast(node);
        }

        for(N successor : _getDirectSuccessors.execute(graph, node))
        {
            processGraph(graph, successor, nodesDFSOrder, storedNodes);
        }

    }
}
