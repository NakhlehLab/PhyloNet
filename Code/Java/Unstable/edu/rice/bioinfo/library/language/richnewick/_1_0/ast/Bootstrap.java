package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/15/11
 * Time: 4:09 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Bootstrap extends AbstractSyntaxNode {

    public <R,T,E extends Exception> R execute(BootstrapAlgo<R,T,E> algo, T input) throws E;
}
