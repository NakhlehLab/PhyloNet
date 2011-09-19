package edu.rice.cs.bioinfo.library.programming;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/18/11
 * Time: 5:40 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Func2<T1,T2,R>
{
    public R execute(T1 input1, T2 input2);
}
