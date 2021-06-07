package edu.rice.cs.bioinfo.library.programming;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/29/13
 * Time: 2:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class Observable3<T1,T2,T3>
{
    private Set<Proc3<T1,T2,T3>> _observers;

    public Observable3()
    {
        _observers = makeSet();
    }

    public void addObserver(Proc3<T1,T2,T3> observer)
    {
        _observers.add(observer);
    }

    public void notify(T1 message1, T2 message2, T3 message3)
    {
        for(Proc3<T1,T2,T3> observer : _observers)
            observer.execute(message1, message2, message3);
    }

    protected Set<Proc3<T1,T2,T3>> makeSet()
    {
        return new HashSet<Proc3<T1,T2,T3>>();
    }
}
