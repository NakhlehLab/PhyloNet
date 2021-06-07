package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/5/13
 * Time: 1:55 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class DiGraphVMBase<N,E>
{
    private final Set<E> _edges;

    public Set<E> getEdges()
    {
        return _edges;
    }

    public DiGraphVMBase(Set<E> edges) {
        _edges = edges;
    }

    public Set<N> getNodes()
    {
        Set<N> nodes = new HashSet<N>();

        for(E edge : _edges)
        {
            nodes.add(getSource(edge));
            nodes.add(getDestination(edge));
        }

        return nodes;
    }

    protected N getEdgeRhs(E edge)
    {
        return getDestination(edge);
    }

    protected N getEdgeLhs(E edge)
    {
        return getSource(edge);
    }

    public abstract N getSource(E edge);

    public abstract N getDestination(E edge);

    public abstract String getNodeLabel(N n);

    public abstract String getEdgeLabel(E edge);


}
