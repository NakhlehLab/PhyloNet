package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/3/11
 * Time: 4:57 PM
 * To change this template use File | Settings | File Templates.
 */
public interface RootageQualifierAlgo<R, T, E extends Exception> {

    public R forEmptyQualifier(RootageQualifierEmpty rootage, T input) throws E;

    public R forNonEmptyQualifier(RootageQualifierNonEmpty rootage, T input) throws E;
}
