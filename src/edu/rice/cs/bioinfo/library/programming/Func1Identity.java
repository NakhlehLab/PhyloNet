package edu.rice.cs.bioinfo.library.programming;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/24/12
 * Time: 2:52 PM
 * To change this template use File | Settings | File Templates.
 */
public final class Func1Identity<T> implements Func1<T,T>
{

    public T execute(T input) {
        return input;
    }
}
