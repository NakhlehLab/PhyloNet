package edu.rice.cs.bioinfo.library.programming;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/18/12
 * Time: 2:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class MutableTuple3<T1,T2,T3> extends MutableTuple<T1,T2>
{
    public T3 Item3;

    public MutableTuple3(T1 item1, T2 item2, T3 item3)
    {
        super(item1, item2);
        Item3 = item3;
    }
}
