package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.jung;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func5;
import edu.uci.ics.jung.graph.DirectedOrderedSparseMultigraph;

import java.math.BigDecimal;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/21/12
 * Time: 3:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphBuilderDirectedOrderedSparse<N,E> extends GraphBuilderJungBase<DirectedOrderedSparseMultigraph<N,E>, N,E>
{
    public GraphBuilderDirectedOrderedSparse(Func1<String,N> makeNode, Func5<N, N, BigDecimal, BigDecimal, BigDecimal, E> makeEdge) {
        super(makeNode, makeEdge);
    }

    @Override
    protected DirectedOrderedSparseMultigraph<N, E> makeEmptyGraph()
    {
        return new DirectedOrderedSparseMultigraph<N, E>();
    }
}