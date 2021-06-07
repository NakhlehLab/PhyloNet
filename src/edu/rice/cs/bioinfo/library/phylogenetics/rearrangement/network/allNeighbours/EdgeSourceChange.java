package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours;

import edu.rice.cs.bioinfo.library.phylogenetics.FindRoot;
import edu.rice.cs.bioinfo.library.phylogenetics.GetDirectPredecessors;
import edu.rice.cs.bioinfo.library.phylogenetics.GetDirectSuccessors;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Tuple;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/11/12
 * Time: 10:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class EdgeSourceChange <G extends Graph<N,E>,N,E> extends NetworkRearrangementOperation<G,N,E>{
            
    public EdgeSourceChange(Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge){
        super(makeNode, makeEdge);
    }

    public G performOperation(){
        Tuple<N,N> nodesOfTargetEdge = _network.getNodesOfEdge(_targetEdge);

        N targetEdgeSibling = null;
        GetDirectSuccessors<N, E> findChildren = new GetDirectSuccessors<N, E>();
        int count = 0;
        for (N node : findChildren.execute(_network, nodesOfTargetEdge.Item1)) {
            if (!node.equals(nodesOfTargetEdge.Item2)) {
                targetEdgeSibling = node;
            }
            count++;
        }
        if(count!=2){
            throw new RuntimeException(nodesOfTargetEdge.Item1+" should have two children!");
        }

        //_originalTargetSibling = targetEdgeSibling;

        N targetEdgeParent = null;
        GetDirectPredecessors<N,E> findParents = new GetDirectPredecessors<N,E>();
        count = 0;
        for(N node : findParents.execute(_network, nodesOfTargetEdge.Item1)){
            targetEdgeParent = node;
            count++;
        }
        if(count!=1 && count!=0){
            throw new RuntimeException(nodesOfTargetEdge.Item1+" should have zero or one parents");
        }

        FindRoot<N> getRoot = new FindRoot<N>();
        N oldRoot = getRoot.execute(_network);

        if (count == 1) {
            if(_network.getEdge(targetEdgeParent, targetEdgeSibling)!=null){
                throw new IllegalStateException();
            }

            _network.removeEdge(_network.getEdge(targetEdgeParent, nodesOfTargetEdge.Item1));
            _network.removeEdge(_network.getEdge(nodesOfTargetEdge.Item1, targetEdgeSibling));
            _sourceEdge = _makeEdge.execute(_network, targetEdgeParent, targetEdgeSibling);
            _network.addEdge(_sourceEdge);

        }
        else{
            _network.removeEdge(_network.getEdge(nodesOfTargetEdge.Item1, targetEdgeSibling));
            _sourceEdge = null;
        }
        
        if(_destinationEdge !=null ){
            Tuple<N,N> nodesOfDestinationEdge = _network.getNodesOfEdge(_destinationEdge);
            E newEdge1OnDestination = _makeEdge.execute(_network, nodesOfDestinationEdge.Item1, nodesOfTargetEdge.Item1);
            E newEdge2OnDestination = _makeEdge.execute(_network, nodesOfTargetEdge.Item1, nodesOfDestinationEdge.Item2);
            _network.addEdge(newEdge1OnDestination);
            _network.addEdge(newEdge2OnDestination);
            _network.removeEdge(_destinationEdge);

        }
        else{
            E newEdge = _makeEdge.execute(_network, nodesOfTargetEdge.Item1, oldRoot);
            _network.addEdge(newEdge);
        }

        return _network;
    }


    public G undoOperation(){
        setParameters(_network, _targetEdge, null, _sourceEdge);
        _network = performOperation();
        return _network;
    }



    /*
    //Not correct for dependent edges
    public void updateNode2Ancestors(Map<N,Set<N>> node2Ancestors){
        Tuple<N,N> nodesOfTargetEdge = _network.getNodesOfEdge(_targetEdge);

        node2Ancestors.get(_originalTargetSibling).remove(nodesOfTargetEdge.Item1);

        Set<N> oldAncestors = node2Ancestors.get(nodesOfTargetEdge.Item1);

        Set<N> newAncestors = new HashSet<N>();
        newAncestors.add(_nodesOfDestinationEdge.Item1);
        newAncestors.addAll(node2Ancestors.get(_nodesOfDestinationEdge.Item1));
        node2Ancestors.put(nodesOfTargetEdge.Item1,newAncestors);

        node2Ancestors.get(_nodesOfDestinationEdge.Item2).add(nodesOfTargetEdge.Item1);
        FindSuccessors<N,E> findAllChildren = new FindSuccessors<N, E>();
        for(N child: findAllChildren.execute(_network, _nodesOfDestinationEdge.Item2)){
            node2Ancestors.get(child).add(nodesOfTargetEdge.Item1);
        }

        Set<N> ancestors = node2Ancestors.get(nodesOfTargetEdge.Item2);
        ancestors.removeAll(oldAncestors);
        ancestors.addAll(newAncestors);
        for(N child: findAllChildren.execute(_network, nodesOfTargetEdge.Item2)){
            ancestors = node2Ancestors.get(child);
            ancestors.removeAll(oldAncestors);
            ancestors.addAll(newAncestors);
        }
    }
    */

}
