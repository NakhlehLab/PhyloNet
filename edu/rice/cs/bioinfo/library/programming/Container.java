package edu.rice.cs.bioinfo.library.programming;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/11
 * Time: 1:24 PM
 * To change this template use File | Settings | File Templates.
 */
public class Container<T> {

    private T _contents;

    public Container(T contents)
    {
        _contents = contents;
    }

    public T getContents()
    {
        return _contents;
    }

    public void setContents(T newValue)
    {
        _contents = newValue;
    }
}
