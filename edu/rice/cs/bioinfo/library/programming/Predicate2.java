package edu.rice.cs.bioinfo.library.programming;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/6/12
 * Time: 4:01 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Predicate2<T1,T2>
{
    public boolean execute(T1 input1, T2 input2);
}
