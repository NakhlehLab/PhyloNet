package edu.rice.cs.bioinfo.programs.soranus.viewModels;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 3/28/13
 * Time: 4:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class VAALOutDataVM implements DocumentVM
{
    public final String Content;

    public VAALOutDataVM(String content) {
        Content = content;
    }

    public <R, E extends Exception> R execute(DocumentVMAlgo<R, E> algo) throws E {
        return algo.forVAALOutDataVM(this);
    }
}