package edu.rice.cs.bioinfo.library.language.richnewick._1_0.graphbuilding.jung;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilder;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.uci.ics.jung.graph.Graph;

import java.io.OptionalDataException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.UUID;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 3:47 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class GraphBuilderJungBase<G extends Graph<N,E>, N,E> implements GraphBuilder<N>
{
    public final G Graph;

    private final Func1<String,N>  _makeNode;

    private final Func5<N,N,BigDecimal,BigDecimal,BigDecimal, E> _makeEdge;

    public GraphBuilderJungBase(Func1<String,N> makeNode, Func5<N,N,BigDecimal,BigDecimal,BigDecimal, E> makeEdge)
    {
        Graph = makeEmptyGraph();
        _makeNode = makeNode;
        _makeEdge = makeEdge;
    }

    protected abstract G makeEmptyGraph();

    public N createNode(String label) {
        N node = _makeNode.execute(label);
        Graph.addVertex(node);
        return node;
    }

    public N createHybridNode(String label, HybridNodeType hybridType, BigInteger hybridNodeIndex) {
        return createNode(label);
    }

    public void createDirectedEdge(N tail, N tip, BigDecimal branchLength, BigDecimal support, BigDecimal probability) {
        E edge = _makeEdge.execute(tail, tip, branchLength, support, probability);
        Graph.addEdge(edge, tail, tip);
    }
}
