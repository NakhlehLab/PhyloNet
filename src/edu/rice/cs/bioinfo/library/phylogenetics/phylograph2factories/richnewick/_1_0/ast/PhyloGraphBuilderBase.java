package edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.GraphBuilder;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;

import java.math.BigInteger;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/6/12
 * Time: 4:14 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class PhyloGraphBuilderBase<G extends Graph<RNNode,?>, D> implements GraphBuilder<RNNode>
{
    public G Graph;

    protected PhyloGraphBuilderBase(G emptyGraph)
    {
        Graph = emptyGraph;
    }

    public RNNode createNode(String label) {

        RNNode node = new RNNode(label, null);
        Graph.addNode(node);
        return node;
    }

    public RNNode createHybridNode(String label, HybridNodeType hybridType, BigInteger hybridNodeIndex) {
        RNNode node = new RNNode(label, hybridType);
        Graph.addNode(node);
        return node;
    }

  //  public abstract void createDirectedEdge(RNNode tail, RNNode tip, BigDecimal branchLength, BigDecimal support, BigDecimal probability);
}
