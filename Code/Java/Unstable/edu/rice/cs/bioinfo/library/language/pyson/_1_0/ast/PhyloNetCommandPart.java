package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 2:20 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class PhyloNetCommandPart extends PySONNodeLineAndCol {

    public PhyloNetCommandPart(int line, int col)
    {
        super(line, col);
    }

    public abstract <R,T,E extends Exception> R execute(PhyloNetCommandPartAlgo<R,T,E> algo, T input) throws E;
}
