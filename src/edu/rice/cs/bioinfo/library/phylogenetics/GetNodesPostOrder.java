package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func1;

import java.util.HashSet;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/10/12
 * Time: 2:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class GetNodesPostOrder<N,E> implements Func1<GraphReadOnly<N,E>, Iterable<N>>{

    private final GetDirectSuccessors<N,E> _getDirectSuccessors = new GetDirectSuccessors<N, E>();

    public Iterable<N> execute(GraphReadOnly<N, E> graph) {

        if(!graph.isRooted())
        {
            throw new IllegalArgumentException("Given graph must be rooted");
        }

        LinkedList<N> nodesPostOrder = new LinkedList<N>();

        processGraph(graph, new FindRoot<N>().execute(graph), nodesPostOrder, new HashSet<N>());

        return nodesPostOrder;
    }

    private void processGraph(GraphReadOnly<N, E> graph, N node, LinkedList<N> nodesPostOrder, HashSet<N> storedNodes) {


        for(N successor : _getDirectSuccessors.execute(graph, node))
        {
            processGraph(graph, successor, nodesPostOrder, storedNodes);
        }

        if(!storedNodes.contains(node))
        {
            nodesPostOrder.addLast(node);
        }
    }
}
