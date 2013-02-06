package edu.rice.cs.bioinfo.programs.soranus.viewModels;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/31/13
 * Time: 6:28 PM
 * To change this template use File | Settings | File Templates.
 */
public interface DocumentVM
{
    public <R,E extends Exception> R execute(DocumentVMAlgo<R,E> algo) throws E;
}
