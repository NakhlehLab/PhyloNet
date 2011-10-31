package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.IdentList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 6:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParameterIdentList extends ParameterBase {

    public Iterable<String> Elements;

    public ParameterIdentList(int line, int column, Iterable<String> elements) {
        super(line, column);

        Elements = elements;
    }

    public <R, T, E extends Exception> R execute(ParameterAlgo<R, T, E> algo, T input) throws E {

        return algo.forIdentList(this, input);
    }
}
