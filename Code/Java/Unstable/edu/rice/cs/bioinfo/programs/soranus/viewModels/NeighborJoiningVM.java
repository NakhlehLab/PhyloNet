package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import javax.swing.*;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/5/13
 * Time: 6:07 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NeighborJoiningVM<N,E> implements DocumentVM
{

    public final Set<E> Edges;

    private N _root = null;

    public N getRoot()
    {
        return _root;
    }

    public void setRoot(N newRoot)
    {
        _root = newRoot;
    }

    public NeighborJoiningVM(Set<E> edges) {
        Edges = edges;
    }

    public abstract Tuple<N,N> getNodesOfEdge(E edge);

    public abstract String getNodeLabel(N node);

    public <R,E extends Exception> R execute(DocumentVMAlgo<R,E> algo) throws E
    {
        return algo.forNeighborJoiningVM(this);
    }


}
