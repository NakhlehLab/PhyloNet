package edu.rice.cs.bioinfo.library.math.discrete;

import edu.rice.cs.bioinfo.library.programming.Proc2;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import javax.management.openmbean.InvalidOpenTypeException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/15/13
 * Time: 1:35 PM
 * To change this template use File | Settings | File Templates.
 */
public class Combinations<T> implements Iterable<Set<T>>
{
    public final int NumToChoose;

    private final T[] _elements;

    private Set<Proc2<Integer,Set<T>>> _elementExhaustedListeneres = new HashSet<Proc2<Integer,Set<T>>>();

    public void addElementExhaustedListener(Proc2<Integer,Set<T>> listener)
    {
        _elementExhaustedListeneres.add(listener);
    }

    public Combinations(Iterable<T> elements, int numToChoose)
    {
        NumToChoose = numToChoose;
        _elements = IterableHelp.toArray(elements);
    }

    public Iterator<Set<T>> iterator()
    {
        final int[] choiceIndexes = new int[NumToChoose];

        for(int i = 0; i<choiceIndexes.length; i++)
        {
            choiceIndexes[i] = i;
        }

        final boolean hasNextInit = NumToChoose <= _elements.length;

        return new Iterator<Set<T>>()
        {
            private boolean _hasNext = hasNextInit;

            private int _lastElementChoiceIndex = choiceIndexes.length-1;

            private boolean _choiceIndexesReflectNewlyExhaustedState = false;

            private int _lastExhaustedIndex = -1;

            public boolean hasNext()
            {
                return _hasNext;
            }

            public Set<T> next()
            {
                Set<T> toBeReturned = new HashSet<T>();
                for(int i = 0; i<choiceIndexes.length; i++)
                {
                    toBeReturned.add(_elements[choiceIndexes[i]]);
                }

                if(_choiceIndexesReflectNewlyExhaustedState)
                {
                    fireElementExhausted(_lastExhaustedIndex, toBeReturned);
                }

                if(choiceIndexes[_lastElementChoiceIndex] != _elements.length -1)
                {
                    choiceIndexes[_lastElementChoiceIndex]++;
                    return toBeReturned;
                }
                else
                {
                    Integer nextAdvancableIndex = findNextFree(_lastElementChoiceIndex);

                    if(nextAdvancableIndex == null)
                    {
                        _hasNext = false;
                        return toBeReturned;
                    }

                    int newValue =  choiceIndexes[nextAdvancableIndex] + 1;
                    choiceIndexes[nextAdvancableIndex] = newValue;
                    int offset = 1;
                    for(int i = nextAdvancableIndex+1; i<choiceIndexes.length; i++)
                    {
                        choiceIndexes[i] = newValue + offset;
                        offset++;
                    }

                    if(nextAdvancableIndex == 0)
                    {
                        _choiceIndexesReflectNewlyExhaustedState = true;
                        _lastExhaustedIndex++;
                    }
                    else
                    {
                        _choiceIndexesReflectNewlyExhaustedState = false;
                    }

                    return toBeReturned;
                }
            }

            private Integer findNextFree(int nonFreeIndex)
            {
                int potentiallyFree = nonFreeIndex -1;
                if(potentiallyFree == -1)
                    return null;

                if(choiceIndexes[potentiallyFree] + 1 == choiceIndexes[nonFreeIndex])
                    return findNextFree(potentiallyFree);
                else
                    return potentiallyFree;
            }

            public void remove()
            {
                throw new InvalidOpenTypeException();
            }
        };
    }

    private void fireElementExhausted(int elementIndex, Set<T> choice)
    {
        for(Proc2<Integer, Set<T>> listener : _elementExhaustedListeneres)
            listener.execute(elementIndex, choice);
    }
}
