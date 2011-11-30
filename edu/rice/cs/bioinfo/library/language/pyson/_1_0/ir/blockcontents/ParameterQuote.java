package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 6:40 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParameterQuote extends ParameterBase {

    public final String UnquotedText;

    public ParameterQuote(int line, int column, String wholeText) {
        super(line, column);
        UnquotedText = wholeText.substring(1, wholeText.length()-1);
    }

    public <R, T, E extends Exception> R execute(ParameterAlgo<R, T, E> algo, T input) throws E {
        return algo.forQuote(this, input);
    }
}
