package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/15/11
 * Time: 4:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class SupportNonEmpty implements Support
{
    public final Text SupportValue;

    public SupportNonEmpty(Text supportValue)
    {
        SupportValue = supportValue;
    }

    public <R, T, E extends Exception> R execute(SupportAlgo<R, T, E> algo, T input) throws E {
        return algo.forSupportNonEmpty(this, input);
    }
}
