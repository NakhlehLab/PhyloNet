package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/3/11
 * Time: 4:55 PM
 * To change this template use File | Settings | File Templates.
 */
public interface RootageQualifier
{
    public <R,T,E extends Exception> R execute(RootageQualifierAlgo<R,T,E> algo, T input) throws E;
}
