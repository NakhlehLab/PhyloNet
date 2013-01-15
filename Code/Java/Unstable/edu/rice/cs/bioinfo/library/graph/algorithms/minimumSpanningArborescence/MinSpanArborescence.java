package edu.rice.cs.bioinfo.library.graph.algorithms.minimumSpanningArborescence;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/15/13
 * Time: 1:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class MinSpanArborescence<E,W>
{
    public final Set<E> Edges;

    public final W SpanWeight;

    public MinSpanArborescence(Set<E> edges, W spanWeight) {
        Edges = edges;
        SpanWeight = spanWeight;
    }
}
