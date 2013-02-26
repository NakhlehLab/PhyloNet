package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge2;
import edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast.NetworkToPhyloGraph2FactoryDefault;
import edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast.RNNode;
import edu.rice.cs.bioinfo.library.programming.Func1Identity;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/6/12
 * Time: 3:57 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("CompareCommonTopologies")
public class CompareCommonTopologies extends CommandBaseFileOut {

    private List<Integer> _unvalidatedNetworksParamPositions = new LinkedList<Integer>();

    private ArrayList<GraphReadOnly<RNNode,PhyloEdge2<RNNode,BigDecimal>>> _validatedInputNetworks = new ArrayList<GraphReadOnly<RNNode,PhyloEdge2<RNNode,BigDecimal>>>(2);

    public CompareCommonTopologies(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                                   Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
        _unvalidatedNetworksParamPositions.add(0);
        _unvalidatedNetworksParamPositions.add(1);
        _validatedInputNetworks.add(null);
        _validatedInputNetworks.add(null);

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

       return checkUnvalidatedNetworkParameters();
    }

    private boolean checkUnvalidatedNetworkParameters()
    {
        boolean noError = true;

        NetworkToPhyloGraph2FactoryDefault graphFactory = new NetworkToPhyloGraph2FactoryDefault(new Func1Identity<BigDecimal>());
        Iterator<Integer> unvalidatedNetworksParamPositions = _unvalidatedNetworksParamPositions.iterator();
        while(unvalidatedNetworksParamPositions.hasNext())
        {
            int paramPosition = unvalidatedNetworksParamPositions.next();
            Network topology =  this.assertAndGetNetwork(paramPosition);
            noError = topology != null && noError;

            if(noError && topology != RuntimeDefinedNetwork.Singleton)
            {
                GraphReadOnly<RNNode,PhyloEdge2<RNNode,BigDecimal>> top = graphFactory.make(topology);
                _validatedInputNetworks.set(paramPosition, top);
                noError = noError && assertAllNodesLabeled(top, (ParameterIdent) this.params.get(paramPosition));
                unvalidatedNetworksParamPositions.remove();
            }
        }

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
        boolean noError = checkUnvalidatedNetworkParameters();

        if(!noError)
            return "";

        BigDecimal largestBLDelta = null;
        BigDecimal largestBLDeltaPercent = null;
        String largestBLDeltaEdge = null;

        BigDecimal largestHybridDelta = null;
        BigDecimal largestHybridDeltaPercent = null;
        String largestHybridDeltaEdge = null;

        StringBuffer result = new StringBuffer();

        GraphReadOnly<RNNode,PhyloEdge2<RNNode,BigDecimal>> topology1 = _validatedInputNetworks.get(0);
         GraphReadOnly<RNNode,PhyloEdge2<RNNode,BigDecimal>> topology2 = _validatedInputNetworks.get(1);

        for(PhyloEdge2<RNNode,BigDecimal> edge1 : topology1.getEdges())
        {

            if(edge1.getProbability() != null)  // skip hybrid edges for jr
                continue;

            PhyloEdge2<RNNode,BigDecimal> edge2 = topology2.getEdge(edge1.Source, edge1.Destination);

            BigDecimal blDelta = edge2.getBranchLength().subtract(edge1.getBranchLength()).abs();

            if(largestBLDelta == null || blDelta.compareTo(largestBLDelta) == 1)
            {
                largestBLDelta = blDelta;
                largestBLDeltaPercent = blDelta.divide(edge1.getBranchLength(), MathContext.DECIMAL128);
                largestBLDeltaEdge = "(" + edge1.Source.Label + ", " + edge1.Destination.Label + ")";
            }

            BigDecimal hybridDelta = null;
            if(edge1.getProbability() != null)
            {

                hybridDelta = edge2.getProbability().subtract(edge1.getProbability()).abs();

                if(largestHybridDelta == null || hybridDelta.compareTo(largestHybridDelta) == 1)
                {
                    largestHybridDelta = hybridDelta;
                    largestHybridDeltaPercent = hybridDelta.divide(edge1.getProbability(), MathContext.DECIMAL128);
                    largestHybridDeltaEdge = "(" + edge1.Source.Label + ", " + edge1.Destination.Label + ")";
                }
            }

             result.append("\n(" + edge1.Source.Label + ", " + edge1.Destination.Label + ")\t" + edge1.getBranchLength() + "\t" + edge2.getBranchLength() +
                           "\t" + blDelta + "\t" + edge1.getProbability() + "\t" + edge2.getProbability() + "\t" + hybridDelta);
        }

        result.append("\nLargest branch length delta: " + largestBLDelta + ": " + largestBLDeltaEdge + ":" + largestBLDeltaPercent + "\n" +
               "Largest hybrid prob   delta: " + largestHybridDelta + ": " + largestHybridDeltaEdge + ":" + largestHybridDeltaPercent +"\n");
        return result.toString();
    }
}
