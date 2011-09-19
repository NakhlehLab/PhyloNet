package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/15/11
 * Time: 4:09 PM
 * To change this template use File | Settings | File Templates.
 */
public final class SupportEmpty implements Support
{
    public static final SupportEmpty Singleton = new SupportEmpty();

    private SupportEmpty()
    {
    }


    public <R, T, E extends Exception> R execute(SupportAlgo<R, T, E> algo, T input) throws E {
        return algo.forSupportEmpty(this, input);
    }
}
