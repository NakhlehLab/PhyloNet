package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 2:56 PM
 * To change this template use File | Settings | File Templates.
 */
public interface HybridNodeQualifierAlgo<R, T, E extends Exception>
{
    public R forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, T input) throws E;

    public R forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, T input) throws E;

    public R forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, T input) throws E;
}
