package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 6:35 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParameterIdent extends ParameterBase {

    public String Content;

    public ParameterIdent(int line, int column, String content) {
        super(line, column);

        Content = content;
    }

    public <R, T, E extends Exception> R execute(ParameterAlgo<R, T, E> algo, T input) throws E {

        return algo.forIdentifier(this, input);
    }
}
