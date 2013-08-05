package edu.rice.cs.bioinfo.library.programming;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/26/13
 * Time: 5:40 PM
 * To change this template use File | Settings | File Templates.
 */
public class Observable1<T>
{
    private Set<Proc1<T>> _observers;

    //public final

    public Observable1()
    {
        _observers = makeSet();
    }

    public void addObserver(Proc1<T> observer)
    {
        _observers.add(observer);
    }

    public void notify(T message)
    {
        for(Proc1<T> observer : _observers)
            observer.execute(message);
    }

    protected Set<Proc1<T>> makeSet()
    {
        return new HashSet<Proc1<T>>();
    }
}
