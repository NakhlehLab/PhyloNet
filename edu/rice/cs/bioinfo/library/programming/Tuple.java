package edu.rice.cs.bioinfo.library.programming;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 4/6/12
 * Time: 6:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class Tuple<T1,T2>
{
    public final T1 Item1;

    public final T2 Item2;

    public Tuple(T1 item1, T2 item2)
    {
        Item1 = item1;
        Item2 = item2;
    }

    @Override
    public boolean equals(Object candidate)
    {
        try
        {
            return equals((Tuple<T1,T2>)candidate);
        }
        catch (ClassCastException e)
        {
            return false;
        }
    }

    public boolean equals(Tuple<T1,T2> candidate)
    {
        if(candidate == null)
            return false;

        return candidate.Item1.equals(Item1) && candidate.Item2.equals(Item2);
    }

    @Override
    public int hashCode()
    {
        return Item1.hashCode() ^ Item2.hashCode();
    }

    public Object other(Object element)
    {
        if(element.equals(Item1))
        {
            return Item2;
        }
        else if(element.equals(Item2))
        {
            return Item1;
        }
        else
        {
            throw new IllegalArgumentException("Given element is not a member of the tuple.");
        }
    }

    @Override
    public String toString()
    {
        return "(" + Item1 + ", " + Item2 + ")";
    }
}
