package edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 5:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class IterableHelp {

    public static Object[] toArray(Iterable elements)
    {
        LinkedList<Object> accum = new LinkedList<Object>();
        for(Object element : elements)
        {
            accum.add(element);
        }

        return accum.toArray();
    }
}
