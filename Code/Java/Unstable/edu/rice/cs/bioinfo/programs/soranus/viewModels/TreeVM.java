package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/5/13
 * Time: 1:58 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TreeVM<N,E> extends GraphVMBase<N,E> implements DocumentVM
{
    public TreeVM(Set<E> edges)
    {
        super(edges);
    }

    public <R,E extends Exception> R execute(DocumentVMAlgo<R,E> algo) throws E
    {
        return algo.forTreeVM(this);
    }
}
