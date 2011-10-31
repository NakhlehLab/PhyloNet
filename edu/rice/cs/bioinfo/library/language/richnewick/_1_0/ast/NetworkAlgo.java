package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/27/11
 * Time: 1:34 PM
 * To change this template use File | Settings | File Templates.
 */
public interface NetworkAlgo<R,T,E extends Exception> {

    public R forNetworkEmpty(NetworkEmpty network, T input) throws E;

    public R forNetworkNonEmpty(NetworkNonEmpty network, T input) throws E;
}
