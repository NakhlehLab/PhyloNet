package edu.rice.cs.bioinfo.library.graph.algorithms.minimumSpanningArborescence;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/14/13
 * Time: 6:39 PM
 * To change this template use File | Settings | File Templates.
 */
public interface MinSpanArborescenceSolver<V,E,W>
{
    public MinSpanArborescence<E,W> tryFindMinSpanArborescence(Set<V> vertices, Set<E> edges);
}
