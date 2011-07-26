package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/15/11
 * Time: 10:49 AM
 * To change this template use File | Settings | File Templates.
 */
public class HybridNodeQualifierWithType extends HybridNodeQualifierNonEmpty
{
    public final Text HybridNodeType;

    public HybridNodeQualifierWithType(Text hybridNodeindex, Text hybridNodeType)
    {
        super(hybridNodeindex);
        HybridNodeType = hybridNodeType;
    }

    @Override
    public <R, T, E extends Exception> R execute(HybridNodeQualifierAlgo<R, T, E> algo, T input) throws E {
        return algo.forHybridNodeQualifierWithType(this, input);
    }
}
