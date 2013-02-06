package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/4/13
 * Time: 7:00 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TransMapVM<N,E> implements DocumentVM
{
    private final Set<E> _edges;

    public Set<E> getEdges()
    {
        return _edges;
    }

    public TransMapVM(Set<E> edges) {
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

    public abstract N getSource(E edge);

    public abstract N getDestination(E edge);

    public abstract String getNodeLabel(N n);

    public <R,E extends Exception> R execute(DocumentVMAlgo<R,E> algo) throws E
    {
        return algo.forTransMapVM(this);
    }


}
