package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.phylolib;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.GraphBuilder;
import edu.rice.cs.bioinfo.library.phylogenetics.DirectedPhyloGraphDefault2;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge2;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloGraph2;
import edu.rice.cs.bioinfo.library.programming.Func1;

import java.math.BigDecimal;
import java.math.BigInteger;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/24/12
 * Time: 2:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphBuilderPhyloGraph2<N,D> implements GraphBuilder<N>
{
    public final PhyloGraph2<N,D> Graph = new DirectedPhyloGraphDefault2<N,D>();

    private final Func1<String,N> _makeNode;

    private final Func1<BigDecimal,D> _makeReal;

    public GraphBuilderPhyloGraph2(Func1<String,N> makeNode, Func1<BigDecimal,D> makeReal)
    {
        _makeNode = makeNode;
        _makeReal = makeReal;
    }

    public N createNode(String label) {
        N node = _makeNode.execute(label);
        Graph.addNode(node);
        return node;
    }

    public N createHybridNode(String label, HybridNodeType hybridType, BigInteger hybridNodeIndex) {
        return createNode(label);
    }

    public void createDirectedEdge(N tail, N tip, BigDecimal branchLength, BigDecimal support, BigDecimal probability) {
        PhyloEdge2<N,D> edge = new PhyloEdge2<N,D>(tail, tip);
        edge.setBranchLength(_makeReal.execute(branchLength));
        edge.setProbability(_makeReal.execute(probability));
        edge.setSupport(_makeReal.execute(support));

        Graph.addEdge(edge);

    }
}
