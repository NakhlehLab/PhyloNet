package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours;

import edu.rice.cs.bioinfo.library.phylogenetics.FindRoot;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Tuple;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/7/12
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReticulationEdgeAddition <G extends Graph<N,E>,N,E> extends NetworkRearrangementOperation<G,N,E>{
    private Tuple<N,N> _nodesOfTargetEdge;

    public ReticulationEdgeAddition(Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge){
        super(makeNode, makeEdge);
    }

    public N getReticulationNode(){
        return _nodesOfTargetEdge.Item2;
    }


    public G performOperation(){
        N sourceEdgeNewNode = _nodesOfTargetEdge==null?_makeNode.execute(_network):_nodesOfTargetEdge.Item1;

        if(_sourceEdge!=null){
            _network.addNode(sourceEdgeNewNode);
            Tuple<N,N> nodesOfSourceEdge = _network.getNodesOfEdge(_sourceEdge);
            _network.removeEdge(_sourceEdge);
            E newEdge1OnSource = _makeEdge.execute(_network, nodesOfSourceEdge.Item1, sourceEdgeNewNode);
            E newEdge2OnSource = _makeEdge.execute(_network, sourceEdgeNewNode, nodesOfSourceEdge.Item2);
            _network.addEdge(newEdge1OnSource);
            _network.addEdge(newEdge2OnSource);
        }else{
            FindRoot<N> getRoot = new FindRoot<N>();
            N oldRoot = getRoot.execute(_network);
            _network.addNode(sourceEdgeNewNode);
            E newEdge = _makeEdge.execute(_network, sourceEdgeNewNode, oldRoot);
            _network.addEdge(newEdge);
        }

        Tuple<N, N> nodesOfDestinationEdge = _network.getNodesOfEdge(_destinationEdge);
        N destinationEdgeNewNode =  _nodesOfTargetEdge==null?_makeNode.execute(_network):_nodesOfTargetEdge.Item2;
        _network.addNode(destinationEdgeNewNode);
        _network.removeEdge(_destinationEdge);
        E newEdge1OnDestination = _makeEdge.execute(_network, nodesOfDestinationEdge.Item1, destinationEdgeNewNode);
        E newEdge2OnDestination = _makeEdge.execute(_network, destinationEdgeNewNode, nodesOfDestinationEdge.Item2);
        _network.addEdge(newEdge1OnDestination);
        _network.addEdge(newEdge2OnDestination);

        _targetEdge = _makeEdge.execute(_network, sourceEdgeNewNode, destinationEdgeNewNode);
        _network.addEdge(_targetEdge);
        _nodesOfTargetEdge = new Tuple<N, N>(sourceEdgeNewNode, destinationEdgeNewNode);
        return _network;
    }


    public G undoOperation(){
        ReticulationEdgeDeletion<G, N, E> undo = new ReticulationEdgeDeletion<G, N, E>(_makeNode, _makeEdge);
        undo.setParameters(_network, _targetEdge, _sourceEdge, _destinationEdge);
        _network = undo.performOperation();
        return _network;
    }


    public void setParameters(G network, E targetEdge, E sourceEdge, E destinationEdge){
        _network = network;
        _targetEdge = targetEdge;
        _sourceEdge = sourceEdge;
        _destinationEdge = destinationEdge;
        _nodesOfTargetEdge = null;
    }

    public void setParameters(G network, E targetEdge, E sourceEdge, E destinationEdge, Tuple<N,N> nodesOfTargetEdge){
        _network = network;
        _targetEdge = targetEdge;
        _sourceEdge = sourceEdge;
        _destinationEdge = destinationEdge;
        _nodesOfTargetEdge = nodesOfTargetEdge;
    }

    /*
    public void updateNode2Ancestors(Map<N,Set<N>> node2Ancestors){
        Tuple<N,N> nodesOfSourceEdge = _network.getNodesOfEdge(_sourceEdge);
        Tuple<N,N> nodesOfDestinationEdge = _network.getNodesOfEdge(_destinationEdge);
        Tuple<N,N> nodesOfTargetEdge = _network.getNodesOfEdge(_targetEdge);
        
        Set<N> targetEdgeSourceAncestors = new HashSet<N>();
        targetEdgeSourceAncestors.addAll(node2Ancestors.get(nodesOfSourceEdge.Item1));
        targetEdgeSourceAncestors.add(nodesOfSourceEdge.Item1);
        node2Ancestors.put(nodesOfTargetEdge.Item1, targetEdgeSourceAncestors);

        FindSuccessors<N,E> findSuccessors = new FindSuccessors<N,E>();
        for(N node : findSuccessors.execute(_network, nodesOfTargetEdge.Item1)){
            node2Ancestors.get(node).add(nodesOfTargetEdge.Item1);
        }

        Set<N> targetEdgeDestinationAncestors = new HashSet<N>();
        targetEdgeDestinationAncestors.addAll(node2Ancestors.get(nodesOfDestinationEdge.Item1));
        targetEdgeDestinationAncestors.add(nodesOfDestinationEdge.Item1);
        targetEdgeDestinationAncestors.addAll(node2Ancestors.get(nodesOfTargetEdge.Item1));
        targetEdgeDestinationAncestors.add(nodesOfTargetEdge.Item1);
        node2Ancestors.put(nodesOfTargetEdge.Item2, targetEdgeDestinationAncestors);

        for(N node : findSuccessors.execute(_network, nodesOfTargetEdge.Item2)){
            Set<N> ancestors = node2Ancestors.get(node);
            ancestors.add(nodesOfTargetEdge.Item2);
            ancestors.add(nodesOfTargetEdge.Item1);
            ancestors.addAll(targetEdgeSourceAncestors);
        }
    }
    */
}
