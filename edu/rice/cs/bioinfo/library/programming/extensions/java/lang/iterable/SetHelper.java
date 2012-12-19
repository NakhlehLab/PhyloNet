package edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable;

import sun.nio.cs.ext.TIS_620;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/6/12
 * Time: 12:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class SetHelper
{
    public static class RelativeComplimentResult<T1,T2>
    {
        public final Set<T1> InAButNotB;

        public final Set<T2> InBButNotA;

        public RelativeComplimentResult(Set<T1> inAButNotB, Set<T2> inBButNotA)
        {
            InAButNotB = inAButNotB;
            InBButNotA = inBButNotA;
        }
    }

    public static  <T1,T2> RelativeComplimentResult<T1,T2> relativeCompliment(Set<T1> a, Set<T2> b)
    {
        RelativeComplimentResult result = new RelativeComplimentResult(new HashSet<T1>(), new HashSet<T2>());

        for(T1 element : a)
        {
            if(!b.contains(a))
            {
                result.InAButNotB.add(a);
            }
        }

        for(T2 element : b)
        {
            if(!a.contains(b))
            {
                result.InBButNotA.add(a);
            }
        }

        return result;
    }
}
