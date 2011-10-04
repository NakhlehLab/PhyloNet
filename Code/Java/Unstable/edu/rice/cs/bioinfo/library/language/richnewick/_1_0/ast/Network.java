package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

import java.nio.ReadOnlyBufferException;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 4:10 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Network extends AbstractSyntaxNode {

    public <R,T,E extends Exception> R execute(NetworkAlgo<R,T,E> algo, T input) throws E;

}
