package edu.rice.cs.bioinfo.library.phylogenetics.phylographfactories.richnewick._1_0.ast;

import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloGraph;
import edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast.PhyloGraphBuilderBase;
import edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast.RNNode;

import java.math.BigDecimal;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/1/12
 * Time: 2:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloGraphBuilder extends PhyloGraphBuilderBase<PhyloGraph<RNNode>,Double>
{
    public PhyloGraph<RNNode> Graph;

    public PhyloGraphBuilder(PhyloGraph<RNNode> emptyGraph)
    {
        super(emptyGraph);
    }

    public void createDirectedEdge(RNNode tail, RNNode tip, BigDecimal branchLength, BigDecimal support, BigDecimal probability) {

        PhyloEdge<RNNode> edge = new PhyloEdge<RNNode>(tail, tip);
        if(branchLength != null)
            edge.setBranchLength(branchLength.doubleValue());
        if(support != null)
            edge.setSupport(support.doubleValue());
        if(probability != null)
            edge.setProbability(probability.doubleValue());
        Graph.addEdge(edge);

    }
}
