package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/7/12
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReticulationFlip extends NetworkRearrangementOperation {

    public boolean performOperation(){
        if(_targetEdge.Item1.isRoot()){
            return false;
        }

        if(_sourceEdge == null){
            _sourceEdge = findParentAndAnotherChild(_targetEdge.Item1, _targetEdge.Item2);
            _destinationEdge = findAnotherParentAndChild(_targetEdge.Item2, _targetEdge.Item1);
            _targetEdgeBrlen = _targetEdge.Item2.getParentDistance(_targetEdge.Item1);
            _targetEdgeInheriProb = _targetEdge.Item2.getParentProbability(_targetEdge.Item1);
        }

        _targetEdge.Item1.removeChild(_targetEdge.Item2);
        _targetEdge.Item2.adoptChild(_targetEdge.Item1, _targetEdgeBrlen);
        _targetEdge = new Tuple<>(_targetEdge.Item2,_targetEdge.Item1);
        _targetEdge.Item2.setParentProbability(_targetEdge.Item1,_targetEdgeInheriProb);
        _targetEdge.Item1.setParentProbability(_destinationEdge.Item1, NetNode.NO_PROBABILITY);
        _targetEdge.Item2.setParentProbability(_sourceEdge.Item1, _targetEdgeInheriProb==NetNode.NO_PROBABILITY? NetNode.NO_PROBABILITY: 1-_targetEdgeInheriProb);

        return true;
    }


    public void undoOperation(){
        ReticulationFlip undo = new ReticulationFlip();
        undo.setParameters(_network, _targetEdge, _targetEdgeBrlen, _targetEdgeInheriProb, _destinationEdge, null, null, _sourceEdge, null, null);
        undo.performOperation();
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
