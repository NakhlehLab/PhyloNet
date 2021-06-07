package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/4/13
 * Time: 7:00 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TransMapVM<N,E> extends DiGraphVMBase<N,E> implements DocumentVM
{

    public TransMapVM(Set<E> edges) {
        super(edges);

    }

    public <R,E extends Exception> R execute(DocumentVMAlgo<R,E> algo) throws E
    {
        return algo.forTransMapVM(this);
    }


}
