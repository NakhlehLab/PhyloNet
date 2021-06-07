package edu.rice.cs.bioinfo.library.programming;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/7/13
 * Time: 5:48 PM
 * To change this template use File | Settings | File Templates.
 */
public interface FuncEx<R,E extends Exception>
{
    public R execute() throws E;
}
