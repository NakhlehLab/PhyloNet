package edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast;

import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge2;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloGraph2;
import edu.rice.cs.bioinfo.library.programming.Func1;

import java.math.BigDecimal;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/1/12
 * Time: 2:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloGraphBuilder<D> extends PhyloGraphBuilderBase<PhyloGraph2<RNNode,D>,D>
{
    private final Func1<BigDecimal,D> _makeD;

    PhyloGraphBuilder(PhyloGraph2<RNNode, D> emptyGraph, Func1<BigDecimal,D> makeD)
    {
        super(emptyGraph);
        _makeD = makeD;
    }

    public void createDirectedEdge(RNNode tail, RNNode tip, BigDecimal branchLength, BigDecimal support, BigDecimal probability) {

        PhyloEdge2<RNNode,D> edge = new PhyloEdge2<RNNode,D>(tail, tip);
        if(branchLength != null)
            edge.setBranchLength(_makeD.execute(branchLength));
        if(support != null)
            edge.setSupport(_makeD.execute(support));
        if(probability != null)
            edge.setProbability(_makeD.execute(probability));
        Graph.addEdge(edge);

    }

}
