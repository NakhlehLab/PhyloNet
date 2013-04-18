package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/5/13
 * Time: 1:55 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class GraphVMBase<N,E>
{
    private final Set<E> _edges;

    public Set<E> getEdges()
    {
        return _edges;
    }

    public GraphVMBase(Set<E> edges) {
        _edges = edges;
    }

    public Set<N> getNodes()
    {
        Set<N> nodes = new HashSet<N>();

        for(E edge : _edges)
        {
            Tuple<N,N> nodesOfEdge = getNodesOfEdge(edge);
            nodes.add(nodesOfEdge.Item1);
            nodes.add(nodesOfEdge.Item2);
        }

        return nodes;
    }

    public abstract Tuple<N,N> getNodesOfEdge(E edge);


    public abstract String getNodeLabel(N n);

    public abstract String getEdgeLabel(E edge);


}
