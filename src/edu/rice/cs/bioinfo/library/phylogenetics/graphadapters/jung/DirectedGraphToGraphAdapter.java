package edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.uci.ics.jung.graph.DirectedGraph;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/6/12
 * Time: 6:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class DirectedGraphToGraphAdapter<N,E> extends JungGraphToGraphAdapterBase<N,E>
{
    public final DirectedGraph<N,E> Graph;

    public DirectedGraphToGraphAdapter(DirectedGraph<N, E> graph, Func1<E, Tuple<N, N>> edgeToTuple) {
        super(graph, edgeToTuple);
        Graph = graph;
    }

    public boolean isRooted() {
        return true;
    }

    public boolean isDestinationNode(N node, E edge) {

        if(Graph.containsEdge(edge))
        {
            return Graph.getDest(edge).equals(node);
        }
        else
        {
            return false;
        }
    }
}
