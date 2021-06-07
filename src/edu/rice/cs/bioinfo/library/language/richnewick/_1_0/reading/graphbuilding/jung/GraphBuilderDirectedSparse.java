package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.jung;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func5;
import edu.uci.ics.jung.graph.DirectedSparseGraph;

import java.math.BigDecimal;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 3:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphBuilderDirectedSparse<N,E> extends GraphBuilderJungBase<DirectedSparseGraph<N,E>, N,E>
{
    public GraphBuilderDirectedSparse(Func1<String,N> makeNode, Func5<N, N, BigDecimal, BigDecimal, BigDecimal, E> makeEdge) {
        super(makeNode, makeEdge);
    }

    @Override
    protected DirectedSparseGraph<N, E> makeEmptyGraph()
    {
        return new DirectedSparseGraph<N, E>();
    }
}
