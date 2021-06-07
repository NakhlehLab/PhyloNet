package edu.rice.cs.bioinfo.library.programming.counters;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/12
 * Time: 4:52 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Counter<T>
{
    public T getCount();

    public void increment();

    public void zero();
}
