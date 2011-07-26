package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/15/11
 * Time: 4:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class BootstrapNonEmpty implements  Bootstrap
{
    public final Text BootstrapValue;

    public BootstrapNonEmpty(Text bootstrapValue)
    {
        BootstrapValue = bootstrapValue;
    }

    public <R, T, E extends Exception> R execute(BootstrapAlgo<R, T, E> algo, T input) throws E {
        return algo.forBootstrapNonEmpty(this, input);
    }
}
