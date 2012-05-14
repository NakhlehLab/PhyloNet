package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func1;

import javax.print.attribute.standard.MediaSize;
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

        processGraph(graph, new FindRoot<N,E>().execute(graph), nodesPostOrder);

        return nodesPostOrder;
    }

    private void processGraph(GraphReadOnly<N, E> graph, N node, LinkedList<N> nodesPostOrder) {


        for(N sucessor : _getDirectSuccessors.execute(graph, node))
        {
            processGraph(graph, sucessor, nodesPostOrder);
        }

        nodesPostOrder.addLast(node);
    }
}
