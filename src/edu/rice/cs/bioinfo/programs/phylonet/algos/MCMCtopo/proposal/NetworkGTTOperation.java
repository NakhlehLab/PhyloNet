package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.Random;

/**
 * Created by dw20 on 06/04/15.
 */
public abstract class NetworkGTTOperation {
    protected Network _network;
    protected Tuple<NetNode,NetNode> _targetEdge;
    protected double _targetEdgeBrlen;
    protected double _targetEdgeInheriProb;
    protected Tuple<NetNode,NetNode> _sourceEdge;
    protected double[] _sourceEdgeBrlens;
    protected double[] _sourceEdgeInheriProbs;
    protected Tuple<NetNode,NetNode> _destinationEdge;
    protected double[] _destinationEdgeBrlens;
    protected double[] _destinationEdgeInheriProbs;
    protected Random _random = new Random();


    public void setParameters(Network network, Tuple<NetNode,NetNode> targetEdge, double targetEdgeBrlen, double targetEdgeInheriProb, Tuple<NetNode,NetNode> sourceEdge, double[] sourceEdgeBrlens, double[] sourceEdgeInheriProbs, Tuple<NetNode,NetNode> destinationEdge, double[] destinationEdgeBrlens, double[] destinationEdgeInheriProbs){
        _network = network;
        _targetEdge = targetEdge;
        _targetEdgeBrlen = targetEdgeBrlen;
        _targetEdgeInheriProb = targetEdgeInheriProb;
        _sourceEdge = sourceEdge;
        _sourceEdgeBrlens = sourceEdgeBrlens;
        _sourceEdgeInheriProbs = sourceEdgeInheriProbs;
        _destinationEdge = destinationEdge;
        _destinationEdgeBrlens = destinationEdgeBrlens;
        _destinationEdgeInheriProbs = destinationEdgeInheriProbs;
    }


    public void setParameters(Network network, Tuple<NetNode,NetNode> targetEdge, Tuple<NetNode,NetNode> sourceEdge, Tuple<NetNode,NetNode> destinationEdge){
        _network = network;
        _targetEdge = targetEdge;
        _targetEdgeBrlen = -1;
        _targetEdgeInheriProb = -1;
        _sourceEdge = sourceEdge;
        _sourceEdgeBrlens = null;
        _sourceEdgeInheriProbs = null;
        _destinationEdge = destinationEdge;
        _destinationEdgeBrlens = null;
        _destinationEdgeInheriProbs = null;
    }


    protected void addNodeToAnEdge(NetNode node, Tuple<NetNode, NetNode> edge, double[] brlens, double[] inheriProbs){
        edge.Item1.adoptChild(node, brlens[0]);
        node.setParentProbability(edge.Item1, inheriProbs[0]);
        edge.Item1.removeChild(edge.Item2);
        node.adoptChild(edge.Item2, brlens[1]);
        edge.Item2.setParentProbability(node, inheriProbs[1]);
        if(node.isNetworkNode() && !Double.isNaN(inheriProbs[0])){
            node.setParentProbability(findAnotherParentAndChild(node, edge.Item1).Item1, 1-inheriProbs[0]);
        }
        if(edge.Item2.isNetworkNode() && !Double.isNaN(inheriProbs[1])){
            edge.Item2.setParentProbability(findAnotherParentAndChild(edge.Item2, node).Item1, 1-inheriProbs[1]);
        }
    }


    protected void randomlyPartitionAnEdge(Tuple<NetNode, NetNode> edge, double[] brlens, double[] inheriProbs){
        double originalBrlen = edge.Item2.getParentDistance(edge.Item1);
        brlens[0] = originalBrlen * Math.random();
        brlens[1] = originalBrlen - brlens[0];
        inheriProbs[0] = Double.NaN;
        inheriProbs[1] = edge.Item2.getParentProbability(edge.Item1);
    }


    protected void removeNodeFromAnEdge(NetNode node, Tuple<NetNode,NetNode> edge, double[] brlens, double[] inheriProbs){
        brlens[1] = edge.Item2.getParentDistance(node);
        inheriProbs[1] = edge.Item2.getParentProbability(node);
        node.removeChild(edge.Item2);
        if(edge.Item1!=null){
            brlens[0] = node.getParentDistance(edge.Item1);
            inheriProbs[0] = node.getParentProbability(edge.Item1);
            edge.Item1.removeChild(node);
            edge.Item1.adoptChild(edge.Item2, brlens[0] + brlens[1]);
            edge.Item2.setParentProbability(edge.Item1, inheriProbs[1]);
        }
    }


    protected Tuple<NetNode,NetNode> findParentAndAnotherChild(NetNode node, NetNode child){
        if(node.getParentCount()>1){
            throw new RuntimeException(node.getName() + " should have zero or one parent");
        }
        if(node.getChildCount()!=2){
            throw new RuntimeException(node.getName() + " should have two children");
        }
        NetNode anotherChild = null;
        for(Object childNode: node.getChildren()){
            if(childNode!=child){
                anotherChild = (NetNode)childNode;
            }
        }
        NetNode parent = null;
        if(node.getParentCount()!=0){
            parent = (NetNode)node.getParents().iterator().next();
        }
        return new Tuple<>(parent, anotherChild);
    }


    protected Tuple<NetNode,NetNode> findAnotherParentAndChild(NetNode node, NetNode parent){
        if(node.getParentCount()!=2){
            throw new RuntimeException(node.getName() + " should have zero or one parent");
        }
        if(node.getChildCount()!=1){
            throw new RuntimeException(node.getName() + " should have one child");
        }
        NetNode anotherParent = null;
        for(Object parentNode: node.getParents()){
            if(parentNode!=parent){
                anotherParent = (NetNode)parentNode;
            }
        }
        NetNode child = (NetNode)node.getChildren().iterator().next();
        return new Tuple<>(anotherParent, child);
    }

    abstract public double getLogHR();

    abstract public boolean performOperation();

    abstract public void undoOperation();
}
