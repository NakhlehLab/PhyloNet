package edu.rice.cs.bioinfo.library.language.dot_2013_1.printing;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/29/13
 * Time: 3:10 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class DigraphDotPrinter<N,E> extends DotPrinterBase<N,E>
{
    public String toDot(Set<N> nodes, Set<E> edges)
    {
        return super.toDot(nodes, edges, "digraph");
    }

    public String toDot(Set<E> edges)
    {
        return super.toDot(edges, "digraph");
    }

    protected String getEdgeString()
    {
        return "->";
    }

    protected N getEdgeRhs(E edge)
    {
        return getDestination(edge);
    }

    protected N getEdgeLhs(E edge)
    {
        return getSource(edge);
    }

    protected abstract N getSource(E edge);

    protected abstract N getDestination(E edge);
}
