package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 4:12 PM
 * To change this template use File | Settings | File Templates.
 */
public interface NodeLabel extends AbstractSyntaxNode
{
    public <R,T,E extends Exception> R execute(NodeLabelAlgo<R,T,E> algo, T input) throws E;
}
