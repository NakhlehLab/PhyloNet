package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/12/13
 * Time: 2:27 PM
 * To change this template use File | Settings | File Templates.
 */
public interface FindAMinimumSpanningTree<G,E,W>
{
    Tuple<Set<E>,W> execute(G graph) throws GraphDisconnectedException;
}
