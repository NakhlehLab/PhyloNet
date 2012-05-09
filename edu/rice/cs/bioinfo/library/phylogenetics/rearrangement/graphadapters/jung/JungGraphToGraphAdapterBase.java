package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.graphadapters.jung;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/6/12
 * Time: 6:15 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class JungGraphToGraphAdapterBase<N,E> implements edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.Graph<N,E>
{
    public final Graph Graph;

    private final Func1<E, Tuple<N,N>> _edgeToTuple;

    public JungGraphToGraphAdapterBase(Graph graph, Func1<E, Tuple<N, N>> edgeToTuple)
    {
        Graph = graph;
        _edgeToTuple = edgeToTuple;
    }

     public void addNode(N node)
    {
        Graph.addVertex(node);
    }

    public void removeNode(N node)
    {
        Graph.removeVertex(node);
    }

    public void addEdge(E edge) {
        Tuple<N,N> nodesOfEdge = _edgeToTuple.execute(edge);
        Graph.addEdge(edge, nodesOfEdge.Item1, nodesOfEdge.Item2);
    }

    public void removeEdge(E edge) {
        Graph.removeEdge(edge);
    }

    public Iterable<N> getNodes() {
        return Graph.getVertices();
    }

    public Iterable<E> getEdges() {
        return Graph.getEdges();
    }

    public Tuple<N, N> getNodesOfEdge(E edge) {
        Pair<N> nodes = Graph.getEndpoints(edge);
        return new Tuple<N, N>(nodes.getFirst(), nodes.getSecond());
    }

      public Iterable<E> getIncidentEdges(N node) {
        return Graph.getIncidentEdges(node);
    }
}
