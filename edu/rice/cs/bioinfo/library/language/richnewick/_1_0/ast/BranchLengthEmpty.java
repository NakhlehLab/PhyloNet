package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/15/11
 * Time: 3:33 PM
 * To change this template use File | Settings | File Templates.
 */
public final class BranchLengthEmpty implements BranchLength{

    public static final BranchLengthEmpty Singleton = new BranchLengthEmpty();

    private BranchLengthEmpty()
    {

    }


    public <R, T, E extends Exception> R execute(BranchLengthAlgo<R, T, E> algo, T input) throws E {
        return algo.forBranchLengthEmpty(this, input);
    }
}
