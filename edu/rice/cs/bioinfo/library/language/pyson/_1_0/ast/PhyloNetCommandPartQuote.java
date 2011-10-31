package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/3/11
 * Time: 5:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloNetCommandPartQuote extends PhyloNetCommandPart
{
    public final String TotalText;

    public PhyloNetCommandPartQuote(String part, int line, int col) {
        super(line, col);
        TotalText = part;
    }

    @Override
    public <R, T, E extends Exception> R execute(PhyloNetCommandPartAlgo<R, T, E> algo, T input) throws E {
       return algo.forQuote(this, input);
    }
}
