package edu.rice.cs.bioinfo.library.language.dot_2013_1.printing;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/5/13
 * Time: 2:22 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class GraphDotPrinter<N,E> extends DotPrinterBase<N,E>
{

    public String toDot(Set<E> edges)
    {
        return super.toDot(edges, "graph");
    }

    protected String getEdgeString()
    {
        return "--";
    }

    protected abstract Tuple<N,N> getNodesOfEdge(E edge);

}
