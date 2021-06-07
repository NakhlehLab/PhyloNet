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
 * Date: 6/7/12
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReticulationEdgeDeletion <G extends Graph<N,E>,N,E> extends NetworkRearrangementOperation<G,N,E>{
    Tuple<N,N> _nodesOfTargetEdge;

    public ReticulationEdgeDeletion(Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge){
        super(makeNode, makeEdge);
    }
    
    public N getReticulationNode(){
        return _nodesOfTargetEdge.Item2;
    }

    public G performOperation(){
        _nodesOfTargetEdge = _network.getNodesOfEdge(_targetEdge);
        N newSourceEdgeItem1 = null, newSourceEdgeItem2 = null, newDestinationEdgeItem1 = null, newDestinationEdgeItem2 = null;
        if(_sourceEdge == null){           
            GetDirectPredecessors<N,E> findParents = new GetDirectPredecessors<N,E>();
            int count = 0;
            for(N node : findParents.execute(_network, _nodesOfTargetEdge.Item1)){
                newSourceEdgeItem1 = node;
                count++;
            }
            if(count!=1 && count!=0){
                throw new RuntimeException(_nodesOfTargetEdge.Item1+" should have zero or one parent");
            }
    
            count = 0;
            for(N node : findParents.execute(_network, _nodesOfTargetEdge.Item2)){
                if(!node.equals(_nodesOfTargetEdge.Item1)){
                    newDestinationEdgeItem1 = node;
                }
                count++;
            }
            if(count!=2){
                throw new RuntimeException(_nodesOfTargetEdge.Item2+" should have two parents");
            }
    
            GetDirectSuccessors<N,E> findChildren = new GetDirectSuccessors<N,E>();
            count = 0;
            for(N node : findChildren.execute(_network, _nodesOfTargetEdge.Item1)){
                if(!node.equals(_nodesOfTargetEdge.Item2)){
                    newSourceEdgeItem2 = node;
                }
                count++;
            }
            if(count!=2){
                throw new RuntimeException(_nodesOfTargetEdge.Item1+" should have two children");
            }
    
            count = 0;
            for(N node : findChildren.execute(_network, _nodesOfTargetEdge.Item2)){
                newDestinationEdgeItem2 = node;
                count++;
            }
            if(count!=1){
                throw new RuntimeException(_nodesOfTargetEdge.Item2+" should have one child");

            }


            if(_network.getEdge(newDestinationEdgeItem1,newDestinationEdgeItem2)!=null){
                throw new IllegalStateException();
            }
            _destinationEdge = _makeEdge.execute(_network, newDestinationEdgeItem1, newDestinationEdgeItem2);

            if(newSourceEdgeItem1!=null){
                _sourceEdge = _makeEdge.execute(_network, newSourceEdgeItem1, newSourceEdgeItem2);
                if(_network.getEdge(newSourceEdgeItem1,newSourceEdgeItem2)!=null || _destinationEdge.equals(_sourceEdge)){
                    throw new IllegalStateException();
                }
                _network.addEdge(_sourceEdge);
            }

            _network.addEdge(_destinationEdge);

        }
        else{
            _network.addEdge(_sourceEdge);
            _network.addEdge(_destinationEdge);
            Tuple<N,N> nodesOfSouceEdge = _network.getNodesOfEdge(_sourceEdge);
            newSourceEdgeItem1 = nodesOfSouceEdge.Item1;
            newSourceEdgeItem2 = nodesOfSouceEdge.Item2;
            Tuple<N,N> nodesOfDestinationEdge = _network.getNodesOfEdge(_destinationEdge);
            newDestinationEdgeItem1 = nodesOfDestinationEdge.Item1;
            newDestinationEdgeItem2 = nodesOfDestinationEdge.Item2;
        }
        //System.out.println("("+newSourceEdgeItem1+","+newSourceEdgeItem2+"),("+newDestinationEdgeItem1+","+newDestinationEdgeItem2+")");
        
        if(newSourceEdgeItem1!=null){
            _network.removeEdge(_network.getEdge(newSourceEdgeItem1, _nodesOfTargetEdge.Item1));
        }

        _network.removeEdge(_network.getEdge(_nodesOfTargetEdge.Item1, newSourceEdgeItem2));
        _network.removeEdge(_network.getEdge(newDestinationEdgeItem1, _nodesOfTargetEdge.Item2));
        _network.removeEdge(_network.getEdge(_nodesOfTargetEdge.Item2, newDestinationEdgeItem2));
        _network.removeEdge(_targetEdge);
        _network.removeNode(_nodesOfTargetEdge.Item1);
        _network.removeNode(_nodesOfTargetEdge.Item2);

        return _network;
    }


    public G undoOperation(){
        ReticulationEdgeAddition<G, N, E> undo = new ReticulationEdgeAddition<G, N, E>(_makeNode, _makeEdge);
        undo.setParameters(_network, null, _sourceEdge, _destinationEdge, _nodesOfTargetEdge);
        _network = undo.performOperation();
        return _network;
    }





    /*
    public void updateNode2Ancestors(Map<N,Set<N>> node2Ancestors){
        Tuple<N,N> nodesOfSourceEdge = _network.getNodesOfEdge(_sourceEdge);
        Tuple<N,N> nodesOfDestinationEdge = _network.getNodesOfEdge(_destinationEdge);


        FindSuccessors<N,E> findSuccessors = new FindSuccessors<N,E>();
        for(N node : findSuccessors.execute(_network, nodesOfSourceEdge.Item1)){
            node2Ancestors.get(node).remove(_nodesOfTargetEdge.Item1);
        }

        for(N node : findSuccessors.execute(_network, nodesOfDestinationEdge.Item1)){
            node2Ancestors.get(node).remove(_nodesOfTargetEdge.Item2);
            node2Ancestors.get(node).remove(_nodesOfTargetEdge.Item1);
            node2Ancestors.get(node).removeAll(node2Ancestors.get(_nodesOfTargetEdge.Item1));
        }

        node2Ancestors.remove(_nodesOfTargetEdge.Item1);
        node2Ancestors.remove(_nodesOfTargetEdge.Item2);
    }
    */
}
