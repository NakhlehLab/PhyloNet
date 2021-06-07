package edu.rice.cs.bioinfo.library.graph.algorithms.bridge;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/10/13
 * Time: 11:30 AM
 * To change this template use File | Settings | File Templates.
 */
public interface IsBridge<G,E>
{
    public boolean execute(G graph, E edge);
}
