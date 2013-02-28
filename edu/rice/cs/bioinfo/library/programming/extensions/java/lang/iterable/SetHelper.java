package edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable;

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
        RelativeComplimentResult<T1,T2> result = new RelativeComplimentResult<T1,T2>(new HashSet<T1>(), new HashSet<T2>());

        for(T1 element : a)
        {
            if(!b.contains(element))
            {
                result.InAButNotB.add(element);
            }
        }

        for(T2 element : b)
        {
            if(!a.contains(element))
            {
                result.InBButNotA.add(element);
            }
        }

        return result;
    }
}
