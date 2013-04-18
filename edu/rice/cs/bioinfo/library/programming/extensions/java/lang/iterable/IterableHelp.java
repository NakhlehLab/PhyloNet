/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable;


import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Predicate1;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 5:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class IterableHelp {

    public static <T> T[] toArray(Iterable<T> elements)
    {
        LinkedList<T> accum = new LinkedList<T>();
        for(T element : elements)
        {
            accum.add(element);
        }
        return (T[])accum.toArray();
    }

    public static <T2,T1 extends T2> List<T2> toList(Iterable<T1> elements)
    {
        List<T2> result = new LinkedList<T2>();

        for(T1 element : elements)
        {
            result.add(element);
        }
        return result;
    }

    public static int countInt(Iterable elements)
    {
        int i = 0;
        for(Object element : elements)
        {
            i++;
        }

        return i;
    }

    public static <T extends Comparable<T>> Collection<T> maxes(Iterable<T> elements)
    {
        LinkedList<T> maxes = new LinkedList<T>();

        for(T element : elements)
        {
            if(maxes.size() == 0)
            {
                maxes.add(element);
            }
            else
            {
                int comparison = maxes.getFirst().compareTo(element);
                if(comparison == 0)
                {
                    maxes.add(element);
                }
                else if(comparison < 0)
                {
                    maxes.clear();
                    maxes.add(element);
                }
            }
        }

        return maxes;
    }

    public static <T extends Comparable<T>> Collection<T> mins(Iterable<T> elements)
    {
        LinkedList<T> mins = new LinkedList<T>();

        for(T element : elements)
        {
            if(mins.size() == 0)
            {
                mins.add(element);
            }
            else
            {
                int comparison = mins.getFirst().compareTo(element);
                if(comparison == 0)
                {
                    mins.add(element);
                }
                else if(comparison > 0)
                {
                    mins.clear();
                    mins.add(element);
                }
            }
        }

        return mins;
    }

    public static <T extends Comparable<T>> T max(Iterable<T> elements)
    {
        return maxes(elements).iterator().next();
    }

    public static <T extends Comparable<T>> T min(Iterable<T> elements)
    {
        return mins(elements).iterator().next();
    }

    public static <T> Collection<T> filter(Iterable<T> elements, Predicate1<T> condition)
    {
        LinkedList<T> filtered = new LinkedList<T>();

        for(T element : elements)
        {
            if(condition.execute(element))
            {
                filtered.add(element);
            }
        }

        return filtered;
    }

    public static Collection<Object> filterUnknown(Iterable<?> elements, Predicate1<Object> condition)
    {
        LinkedList<Object> filtered = new LinkedList<Object>();

        for(Object element : elements)
        {
            if(condition.execute(element))
            {
                filtered.add(element);
            }
        }

        return filtered;
    }

    public static <T extends Comparable> Iterable<T> sortAscending(Iterable<T> elements)
    {
        List<T> elementsList = toList(elements);
        Collections.sort(elementsList);
        return elementsList;
    }

    public static <T,R> Iterable<R> map(final Iterable<T> elements, final Func1<T,R> func)
    {
        return new Iterable<R>()
        {
            public Iterator<R> iterator()
            {
                return new Iterator<R>()
                {
                    private Iterator<T> _elements = elements.iterator();

                    public boolean hasNext()
                    {
                        return _elements.hasNext();
                    }

                    public R next()
                    {
                        return func.execute(_elements.next());
                    }

                    public void remove()
                    {
                        _elements.remove();
                    }
                };
            }
        };
    }
}
