package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/3/11
 * Time: 5:42 PM
 * To change this template use File | Settings | File Templates.
 */
public class RootageQualifierNonEmpty implements RootageQualifier
{
    public String Qualifier;

    public RootageQualifierNonEmpty(String qualifier)
    {
        Qualifier = qualifier;
    }

    public <R, T, E extends Exception> R execute(RootageQualifierAlgo<R, T, E> algo, T input) throws E {
        return algo.forNonEmptyQualifier(this, input);
    }
}
