package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/15/11
 * Time: 10:55 AM
 * To change this template use File | Settings | File Templates.
 */
public class HybridNodeQualifierNonEmpty implements HybridNodeQualifier
{
    public final Text HybridNodeIndex;

    public HybridNodeQualifierNonEmpty(Text hybridNodeIndex)
    {
       HybridNodeIndex = hybridNodeIndex;
    }

    public <R, T, E extends Exception> R execute(HybridNodeQualifierAlgo<R, T, E> algo, T input) throws E {
        return algo.forHybridNodeQualifierNonEmpty(this, input);
    }
}
