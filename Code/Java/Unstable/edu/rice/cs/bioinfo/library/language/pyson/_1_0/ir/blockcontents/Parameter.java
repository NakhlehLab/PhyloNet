package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:33 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Parameter {

    public int getLine();

    public int getColumn();

    public <R,T,E extends Exception> R execute(ParameterAlgo<R,T,E> algo, T input) throws E;
}
