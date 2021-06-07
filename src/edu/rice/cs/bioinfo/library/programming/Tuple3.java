package edu.rice.cs.bioinfo.library.programming;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/18/12
 * Time: 2:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class Tuple3<T1,T2,T3> extends Tuple<T1,T2>
{
    public final T3 Item3;

    public Tuple3(T1 item1, T2 item2, T3 item3)
    {
        super(item1, item2);
        Item3 = item3;
    }

    @Override
    public boolean equals(Object candidate)
    {
        try
        {
            return equals((Tuple3<T1,T2,T3>)candidate);
        }
        catch (ClassCastException e)
        {
            return false;
        }
    }

    public boolean equals(Tuple3<T1,T2,T3> candidate)
    {
        if(candidate == null)
            return false;

        return candidate.Item1.equals(Item1) && candidate.Item2.equals(Item2) && candidate.Item3.equals(Item3);
    }

    @Override
    public int hashCode()
    {
        return Item1.hashCode() ^ Item2.hashCode() ^ Item3.hashCode();
    }

    @Override
    public String toString()
    {
        return "(" + Item1 + ", " + Item2 + ", " + Item3 + ")";
    }
}
