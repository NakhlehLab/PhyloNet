package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.phylogenetics.AreSameTopology;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge2;
import edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast.NetworkToPhyloGraph2FactoryDefault;
import edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast.RNNode;
import edu.rice.cs.bioinfo.library.programming.Func1Identity;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.SetHelper;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/6/12
 * Time: 3:57 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("BranchLengthHybridProbAccumDiff")
public class BranchLengthHybridProbAccumDiff extends CommandBaseFileOut {

    private GraphReadOnly<RNNode,PhyloEdge2<RNNode,BigDecimal>> _network1 = null;

    private GraphReadOnly<RNNode,PhyloEdge2<RNNode,BigDecimal>> _network2 = null;

    private Map<RNNode,RNNode> _network1NodeToNetwork2Node;

    public BranchLengthHybridProbAccumDiff(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                                           Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);

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

       NetworkToPhyloGraph2FactoryDefault graphFactory = new NetworkToPhyloGraph2FactoryDefault(new Func1Identity<BigDecimal>());

       boolean noError = true;
       if(_network1 == null)
       {
           Network network1 = this.assertAndGetNetwork(0);
           noError = noError && network1 != null;

           if(network1 != null && network1 != RuntimeDefinedNetwork.Singleton)
           {
               _network1 = graphFactory.make(network1);
               noError = noError && assertAllNodesLabeled(_network2, this.assertParameterIdent(0));
           }
       }

       if(_network2 == null)
       {
           Network network2 = this.assertAndGetNetwork(1);
           noError = noError && network2 != null;

           if(network2 != null && network2 != RuntimeDefinedNetwork.Singleton)
           {
               _network2 = graphFactory.make(network2);
               noError = noError && assertAllNodesLabeled(_network2, this.assertParameterIdent(1));
           }
       }

       if(_network1 != null && _network2 != null)
       {
           int network1EdgesCount = IterableHelp.countInt(_network1.getEdges());
           int network2EdgesCount = IterableHelp.countInt(_network2.getEdges());

           if(network1EdgesCount != network2EdgesCount)
           {
               this.errorDetected.execute("Given networks must have the same number of edges.", this.params.get(0).getLine(), this.params.get(0).getColumn());
               noError = false;
           }

           Map<String,RNNode> labelToNode1 = makeLabelToNode(_network1);
           Map<String,RNNode> labelToNode2 = makeLabelToNode(_network2);

           SetHelper.RelativeComplimentResult<String,String> compliments = SetHelper.relativeCompliment(labelToNode1.keySet(), labelToNode2.keySet());

           if(compliments.InAButNotB.size() == 0 && compliments.InBButNotA.size() == 0)
           {
               Map<RNNode,RNNode> _network1NodeToNetwork2Node = new HashMap<RNNode,RNNode>();

               for(Map.Entry<String,RNNode> labelToNode : labelToNode1.entrySet())
               {
                   _network1NodeToNetwork2Node.put(labelToNode.getValue(), labelToNode2.get(labelToNode.getKey()));
               }


               AreSameTopology.Result sameTop = new AreSameTopology().execute(_network1, _network2, _network1NodeToNetwork2Node);

               if(!sameTop.SameTopology)
               {
                   noError = false;
                   this.errorDetected.execute("Given networks are not of the same topology.", this.getDefiningSyntaxCommand().getLine(), this.getDefiningSyntaxCommand().getColumn());
               }
           }
           else
           {
                noError = false;
               this.errorDetected.execute("Network labels are not one to one. e.g. " + (compliments.InAButNotB.size() > 0 ?
                       compliments.InAButNotB.iterator().next() : compliments.InBButNotA.iterator().next()), this.getDefiningSyntaxCommand().getLine(),
                                                                                                             this.getDefiningSyntaxCommand().getColumn());
           }

       }

       return noError;
    }

    private Map<String, RNNode> makeLabelToNode(GraphReadOnly<RNNode, PhyloEdge2<RNNode, BigDecimal>> network) {

        Map<String,RNNode> labelToNode = new HashMap<String, RNNode>();

        for(RNNode node : network.getNodes())
        {
            labelToNode.put(node.Label, node);
        }

        return labelToNode;

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
        boolean noError = checkParamsForCommand();

        if(!noError)
            return "";

        BigDecimal diffAccum = BigDecimal.ZERO;

        for(PhyloEdge2<RNNode,BigDecimal> edge : _network1.getEdges())
        {
            Tuple<RNNode,RNNode> edgeNodes = _network1.getNodesOfEdge(edge);
            PhyloEdge2<RNNode,BigDecimal> correspondingEdge = _network2.getEdge(_network1NodeToNetwork2Node.get(edgeNodes.Item1),
                                                                                _network1NodeToNetwork2Node.get(edgeNodes.Item2));

            BigDecimal blDelta = edge.getBranchLength().subtract(correspondingEdge.getBranchLength()).abs();
            diffAccum = diffAccum.add(blDelta);

            BigDecimal hybridProbDelta = edge.getProbability().subtract(correspondingEdge.getProbability()).abs();
            diffAccum = diffAccum.add(hybridProbDelta);

        }

        return diffAccum.toString();

    }
}
