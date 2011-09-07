package edu.rice.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 5:44 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Block extends PySONNode{

    public <R,T,E extends Exception> R execute(BlockAlgo<R,T,E> algo, T input) throws E;
}
