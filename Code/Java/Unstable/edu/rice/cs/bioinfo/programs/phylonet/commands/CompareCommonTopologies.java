package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge2;
import edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast.NetworkToPhyloGraph2FactoryDefault;
import edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast.RNNode;
import edu.rice.cs.bioinfo.library.programming.Func1Identity;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/6/12
 * Time: 3:57 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("comparecommontopologies")
public class CompareCommonTopologies extends CommandBaseFileOut {

     private GraphReadOnly<RNNode,PhyloEdge2<RNNode,BigDecimal>> _topology1;

    private GraphReadOnly<RNNode,PhyloEdge2<RNNode,BigDecimal>> _topology2;

    CompareCommonTopologies(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 3;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean  noError = true;

        Network topology1 =  this.assertAndGetNetwork(0);
        noError = topology1 != null && noError;

        Network topology2 = this.assertAndGetNetwork(1);
        noError = topology2 != null && noError;

        if(topology1 != null && topology2 != null)
        {
            NetworkToPhyloGraph2FactoryDefault graphFactory = new NetworkToPhyloGraph2FactoryDefault(new Func1Identity<BigDecimal>());
           _topology1 = graphFactory.make(topology1);
           _topology2 = graphFactory.make(topology2);
        }

        noError = noError && assertAllNodesLabeled(_topology1, (ParameterIdent) this.params.get(0));
        noError = noError && assertAllNodesLabeled(_topology2, (ParameterIdent) this.params.get(1));

        return noError;
    }

    private boolean assertAllNodesLabeled(GraphReadOnly<RNNode, PhyloEdge2<RNNode, BigDecimal>> graph, ParameterIdent parameter) {

        for(RNNode node : graph.getNodes())
        {
            if(node.Label == null)
            {
                this.errorDetected.execute("All nodes of network " + parameter.Content + " must be labeled.", parameter.getLine(), parameter.getColumn());
                return false;
            }
        }

        return true;
    }

    @Override
    protected String produceResult()
    {
        BigDecimal largestBLDelta = null;
        String largestBLDeltaEdge = null;

        BigDecimal largestHybridDelta = null;
        String largestHybridDeltaEdge = null;

        StringBuffer result = new StringBuffer();

        for(PhyloEdge2<RNNode,BigDecimal> edge1 : _topology1.getEdges())
        {

            PhyloEdge2<RNNode,BigDecimal> edge2 = _topology2.getEdge(edge1.Source, edge1.Destination);

            BigDecimal blDelta = edge2.getBranchLength().subtract(edge1.getBranchLength()).abs();

            if(largestBLDelta == null || blDelta.compareTo(largestBLDelta) == 1)
            {
                largestBLDelta = blDelta;
                largestBLDeltaEdge = "(" + edge1.Source.Label + ", " + edge1.Destination.Label + ")";
            }

            BigDecimal hybridDelta = null;
            if(edge1.getProbability() != null)
            {

                hybridDelta = edge2.getProbability().subtract(edge1.getProbability()).abs();

                if(largestHybridDelta == null || hybridDelta.compareTo(largestHybridDelta) == 1)
                {
                    largestHybridDelta = hybridDelta;
                    largestHybridDeltaEdge = "(" + edge1.Source.Label + ", " + edge1.Destination.Label + ")";
                }
            }

             result.append("\n(" + edge1.Source.Label + ", " + edge1.Destination.Label + ")\t" + edge1.getBranchLength() + "\t" + edge2.getBranchLength() +
                           "\t" + blDelta + "\t" + edge1.getProbability() + "\t" + edge2.getProbability() + "\t" + hybridDelta);
        }

        result.append("\nLargest branch length delta: " + largestBLDelta + ": " + largestBLDeltaEdge + "\n" +
               "Largest hybrid prob   delta: " + largestHybridDelta + ": " + largestHybridDeltaEdge + "\n");
        return result.toString();
    }
}
