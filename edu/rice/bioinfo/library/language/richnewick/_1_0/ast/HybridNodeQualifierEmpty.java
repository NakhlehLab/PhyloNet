package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/15/11
 * Time: 10:56 AM
 * To change this template use File | Settings | File Templates.
 */
public final class HybridNodeQualifierEmpty implements HybridNodeQualifier
{
    public static final HybridNodeQualifierEmpty Singleton = new HybridNodeQualifierEmpty();

    private HybridNodeQualifierEmpty()
    {
    }


    public <R, T, E extends Exception> R execute(HybridNodeQualifierAlgo<R, T, E> algo, T input) throws E {
        return algo.forHybridNodeQualifierEmpty(this, input);
    }
}
