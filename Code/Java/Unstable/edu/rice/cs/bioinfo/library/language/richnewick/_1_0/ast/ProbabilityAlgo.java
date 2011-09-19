package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 3:13 PM
 * To change this template use File | Settings | File Templates.
 */
public interface ProbabilityAlgo<R, T, E extends Exception>
{
    public R forProbabilityEmpty(ProbabilityEmpty prob, T input) throws E;

    public R forProbabilityNonEmpty(ProbabilityNonEmpty prob, T input) throws E;

}
