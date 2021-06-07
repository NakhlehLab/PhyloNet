package edu.rice.cs.bioinfo.library.programming.counters;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/12
 * Time: 4:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class CounterInt implements Counter<Integer>
{
    private int _count = 0;

    public Integer getCount() {
        return _count;
    }

    public void zero()
    {
        _count = 0;
    }

    public void increment() {
       _count++;
    }
}
