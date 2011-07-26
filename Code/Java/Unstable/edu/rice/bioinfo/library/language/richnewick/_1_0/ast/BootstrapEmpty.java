package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/15/11
 * Time: 4:09 PM
 * To change this template use File | Settings | File Templates.
 */
public final class BootstrapEmpty implements Bootstrap
{
    public static final BootstrapEmpty Singleton = new BootstrapEmpty();

    private BootstrapEmpty()
    {
    }


    public <R, T, E extends Exception> R execute(BootstrapAlgo<R, T, E> algo, T input) throws E {
        return algo.forBootstrapEmpty(this, input);
    }
}
