package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours;

import edu.rice.cs.bioinfo.library.phylogenetics.GetDirectPredecessors;
import edu.rice.cs.bioinfo.library.phylogenetics.GetDirectSuccessors;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Tuple;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/12/12
 * Time: 4:30 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReticulationEdgeDestinationChange<G extends Graph<N,E>,N,E> extends NetworkRearrangementOperation<G,N,E>{

    public ReticulationEdgeDestinationChange(Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge){
        super(makeNode, makeEdge);
    }

    public G performOperation(){
        Tuple<N,N> nodesOfTargetEdge = _network.getNodesOfEdge(_targetEdge);

        N targetEdgeChild = null;
        GetDirectSuccessors<N, E> findChildren = new GetDirectSuccessors<N, E>();
        int count = 0;
        for (N node : findChildren.execute(_network, nodesOfTargetEdge.Item2)) {
            targetEdgeChild = node;
            count++;
        }
        if(count!=1){
            throw new RuntimeException(nodesOfTargetEdge.Item2+" should have one child, not " + count);
        }

        N targetEdgeSibling = null;
        GetDirectPredecessors<N,E> findParents = new GetDirectPredecessors<N,E>();
        count = 0;
        for(N node : findParents.execute(_network, nodesOfTargetEdge.Item2)){
            if(!node.equals(nodesOfTargetEdge.Item1)){
                targetEdgeSibling = node;
            }
            count++;
        }
        if(count!=2){
            throw new RuntimeException(nodesOfTargetEdge.Item2+" should have two parents");
        }

        if(_sourceEdge == null){
            if(_network.getEdge(targetEdgeSibling,targetEdgeChild)!=null){
                throw new IllegalStateException();
            }
            _sourceEdge = _makeEdge.execute(_network, targetEdgeSibling, targetEdgeChild);
        }
        _network.addEdge(_sourceEdge);
        E removeEdge1OnSource = _network.getEdge(targetEdgeSibling, nodesOfTargetEdge.Item2);
        _network.removeEdge(removeEdge1OnSource);
        E removeEdge2OnSource = _network.getEdge(nodesOfTargetEdge.Item2, targetEdgeChild);
         _network.removeEdge(removeEdge2OnSource);


        Tuple<N,N> nodesOfDestinationEdge = _network.getNodesOfEdge(_destinationEdge);
        E newEdge1OnDestination = _makeEdge.execute(_network, nodesOfDestinationEdge.Item1, nodesOfTargetEdge.Item2);
        E newEdge2OnDestination = _makeEdge.execute(_network, nodesOfTargetEdge.Item2, nodesOfDestinationEdge.Item2);
        _network.removeEdge(_destinationEdge);
        _network.addEdge(newEdge1OnDestination);
        _network.addEdge(newEdge2OnDestination);


        return _network;
    }

    public G undoOperation(){
        E temp = _sourceEdge;
        _sourceEdge = _destinationEdge;
        _destinationEdge = temp;
        _network = performOperation();

        return _network;
    }
    



}
