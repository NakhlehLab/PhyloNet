package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 4:18 PM
 * To change this template use File | Settings | File Templates.
 */
public interface BootstrapAlgo<R, T, E extends Exception> {

    public R forBootstrapNonEmpty(BootstrapNonEmpty bootstrap, T input) throws E;

    public R forBootstrapEmpty(BootstrapEmpty bootstrap, T input) throws E;
}
