package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal;

/**
 * Move the head of a random edge from current network
 * Created by dw20 on 5/31/15.
 * parameters: _targetEdge _destinationEdge
 */
public class MoveHead extends NetworkGTTOperation{

    private double hastingRatio = Double.MIN_VALUE;

    public double getLogHR() {
        return hastingRatio;
    }

    public boolean performOperation(){
        if(_sourceEdge == null){
            _sourceEdge = findAnotherParentAndChild(_targetEdge.Item2, _targetEdge.Item1);
            if(_sourceEdge.Item2.hasParent(_sourceEdge.Item1)){
                return false;
            }
        }

        _sourceEdgeBrlens = new double[2];
        _sourceEdgeInheriProbs = new double[2];
        removeNodeFromAnEdge(_targetEdge.Item2, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);

        // dest / source
        hastingRatio = Math.log( _destinationEdge.Item2.getParentDistance(_destinationEdge.Item1) /
                _sourceEdge.Item2.getParentDistance(_sourceEdge.Item1) );

        if(_destinationEdgeBrlens==null){
            _destinationEdgeBrlens = new double[2];
            _destinationEdgeInheriProbs = new double[2];
            randomlyPartitionAnEdge(_destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
        }
        addNodeToAnEdge(_targetEdge.Item2, _destinationEdge,_destinationEdgeBrlens, _destinationEdgeInheriProbs);
        _targetEdge.Item2.setParentProbability(_destinationEdge.Item1, 1-_targetEdge.Item2.getParentProbability(_targetEdge.Item1));

        return true;
    }

    public void undoOperation(){
        MoveHead undo = new MoveHead();
        undo.setParameters(_network, _targetEdge, -1, -1, _destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);
        undo.performOperation();
    }
    



}
