package edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.uci.ics.jung.graph.Graph;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/6/12
 * Time: 5:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class UndirectedGraphToGraphAdapter<N,E> extends JungGraphToGraphAdapterBase<N,E>
{

    public UndirectedGraphToGraphAdapter(Graph<N, E> graph, Func1<E, Tuple<N, N>> edgeToTuple) {
        super(graph, edgeToTuple);
    }

    public boolean isRooted() {
        return false;
    }

    public boolean isDestinationNode(N node, E edge) {
        return false;
    }
}
