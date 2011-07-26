package edu.rice.bioinfo.library.programming;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/25/11
 * Time: 5:34 PM
 * To change this template use File | Settings | File Templates.
 */
public final class Func1Null<T,R> implements  Func1<T,R>
{
    public static final Func1Null Singleton = new Func1Null();

    public R execute(T input) {
        return null;
    }
}
