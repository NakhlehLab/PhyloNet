package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 4:12 PM
 * To change this template use File | Settings | File Templates.
 */
public interface BranchLengthAlgo<R, T, E extends Exception>
{
    public R forBranchLengthEmpty(BranchLengthEmpty branchLength, T input) throws E;

    public R forBranchLengthNonEmpty(BranchLengthNonEmpty branchLength, T input) throws E;
}
