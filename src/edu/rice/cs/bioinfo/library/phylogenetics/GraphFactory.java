package edu.rice.cs.bioinfo.library.phylogenetics;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/1/12
 * Time: 1:40 PM
 * To change this template use File | Settings | File Templates.
 */
public interface GraphFactory<N,E,T>
{
    public Graph<N,E> make(T input);
}
