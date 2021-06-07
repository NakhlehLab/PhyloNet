package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;

import java.util.Random;

/**
 * Created by Yun Yu
 *
 * This class is for network rearrangement
 */
public abstract class NetworkRearrangementOperation {
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


    /**
     * This function is set the parameters
     *
     * @param network               the species network to be rearranged
     * @param targetEdge            the target edge that is involved in this rearrangement
     * @param targetEdgeBrlen       the branch length of <code>targetEdge</code>
     * @param targetEdgeInheriProb  the inheritance probability of <code>targetEdge</code>
     * @param sourceEdge            the source edge that is involved in this rearrangement
     * @param sourceEdgeBrlens      the branch lengths of <code>sourceEdge</code>
     * @param sourceEdgeInheriProbs the inheritance probabilities of <code>sourceEdge</code>
     * @param destinationEdge       the source edge that is involved in this rearrangement
     * @param destinationEdgeBrlens      the branch lengths of <code>destinationEdge</code>
     * @param destinationEdgeInheriProbs the inheritance probabilities of <code>destinationEdge</code>
     */
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


    /**
     * This function is set the parameters
     *
     * @param network               the species network to be rearranged
     * @param targetEdge            the target edge that is involved in this rearrangement
     * @param sourceEdge            the source edge that is involved in this rearrangement
     * @param destinationEdge       the source edge that is involved in this rearrangement
     */
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


    /**
     * This function is to insert a node to an edge
     *
     * @param node          the node to be inserted
     * @param edge          the edge that the node is to be inserted on
     * @param brlens        the branch lengths of the two edges after insertion
     * @param inheriProbs   the inheritance probabilities of the two edges after insertion
     */
    protected void addNodeToAnEdge(NetNode node, Tuple<NetNode, NetNode> edge, double[] brlens, double[] inheriProbs){
        edge.Item1.adoptChild(node, brlens[0]);
        node.setParentProbability(edge.Item1, inheriProbs[0]);
        edge.Item1.removeChild(edge.Item2);
        node.adoptChild(edge.Item2, brlens[1]);
        edge.Item2.setParentProbability(node, inheriProbs[1]);
        if(node.isNetworkNode() && inheriProbs[0]!=NetNode.NO_PROBABILITY){
            node.setParentProbability(findAnotherParentAndChild(node, edge.Item1).Item1, 1-inheriProbs[0]);
        }
        if(edge.Item2.isNetworkNode() && inheriProbs[1]!=NetNode.NO_PROBABILITY){
            edge.Item2.setParentProbability(findAnotherParentAndChild(edge.Item2, node).Item1, 1-inheriProbs[1]);
        }
    }


    /**
     * This function is to randomly partition an edge
     *
     * @param edge          the edge that is to be randomly partitioned
     * @param brlens        the resulting branch lengths after partitioning
     * @param inheriProbs   the resulting inheritance probabilities of the two edges after partitioning
     */
    protected void randomlyPartitionAnEdge(Tuple<NetNode, NetNode> edge, double[] brlens, double[] inheriProbs){
        double originalBrlen = edge.Item2.getParentDistance(edge.Item1);
        if(originalBrlen == NetNode.NO_DISTANCE){
            brlens[0] = NetNode.NO_DISTANCE;
            brlens[1] = NetNode.NO_DISTANCE;
        }
        else {
            brlens[0] = originalBrlen * Math.random();
            brlens[1] = originalBrlen - brlens[0];
        }
        inheriProbs[0] = NetNode.NO_PROBABILITY;
        inheriProbs[1] = edge.Item2.getParentProbability(edge.Item1);
    }


    /**
     * This function is to remove a node from an edge
     *
     * @param node          the node to be removed
     * @param edge          the edge that the node is to be removed from
     * @param brlens        the original branch lengths before removing
     * @param inheriProbs   the original inheritance probabilities before removing
     */
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


    /**
     * This function is to find the parent node and another child node of a given node
     *
     * @param node      the given node
     * @param child     one of the children of <code>node</code>
     *
     * @return  a tuple which contains the parent of <code>node</code> and the child of <code>node</code> which is not <code>child</code>
     */
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


    /**
     * This function is to find the child node and another parent node of a given node
     *
     * @param node      the given node
     * @param parent    one of the parents of <code>node</code>
     *
     * @return  a tuple which contains the parent of <code>node</code> which is not <code>parent</code> and the child of <code>node</code>
     */
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


    /**
     * This function is to perform the operation
     *
     * @return  whether the operation is performed successfully
     *          returns false when the resulting network is invalid
     */
    abstract public boolean performOperation();


    /**
     * This function is to undo the operation
     */
    abstract public void undoOperation();

    public Tuple<NetNode, NetNode> getTargetEdge(){
        return _targetEdge;
    }
    public Tuple<NetNode, NetNode> getSourceEdge(){
        return _sourceEdge;
    }
    public Tuple<NetNode, NetNode> getDestiniationEdge(){
        return _destinationEdge;
    }
}
