package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/15/11
 * Time: 10:43 AM
 * To change this template use File | Settings | File Templates.
 */
public interface HybridNodeQualifier extends AbstractSyntaxNode
{
    public <R,T,E extends Exception> R execute(HybridNodeQualifierAlgo<R,T,E> algo, T input) throws E;
}
