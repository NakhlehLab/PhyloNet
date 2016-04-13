package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * File a reticulation edge from current network
 * Created by dw20 on 5/31/15.
 * parameters: _targetEdge _destEdge
 */
public class FlipReticulation extends NetworkGTTOperation {

    private double hastingRatio = Double.MIN_VALUE;

    // source / dest
    public double getLogHR() {
        return hastingRatio;
    }

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
        double srcDist  = _sourceEdge.Item2.getParentDistance(_targetEdge.Item1) + _targetEdge.Item1.getParentDistance(_sourceEdge.Item1);
        double destDist = _destinationEdge.Item2.getParentDistance(_targetEdge.Item2) +
                _targetEdge.Item2.getParentDistance(_destinationEdge.Item1);
        double newSrcDist = srcDist * _random.nextDouble();
        double newDestDist = destDist * _random.nextDouble();
        hastingRatio = 0.0; //Math.log(srcDist / destDist);

        _targetEdge.Item1.removeChild(_targetEdge.Item2);
        _targetEdge.Item2.adoptChild(_targetEdge.Item1, _targetEdgeBrlen);
        _targetEdge = new Tuple<>(_targetEdge.Item2,_targetEdge.Item1);
        _targetEdge.Item2.setParentProbability(_targetEdge.Item1,_targetEdgeInheriProb);
        _targetEdge.Item1.setParentProbability(_destinationEdge.Item1, NetNode.NO_PROBABILITY);
        _targetEdge.Item2.setParentProbability(_sourceEdge.Item1, 1-_targetEdgeInheriProb);

        _targetEdge.Item2.setParentDistance(_sourceEdge.Item1, newSrcDist);
        _sourceEdge.Item2.setParentDistance(_targetEdge.Item2, srcDist - newSrcDist);

        _targetEdge.Item1.setParentDistance(_destinationEdge.Item1, newDestDist);
        _destinationEdge.Item2.setParentDistance(_targetEdge.Item1, destDist - newDestDist);
        return true;
    }


    public void undoOperation(){
        FlipReticulation undo = new FlipReticulation();
        undo.setParameters(_network, _targetEdge, _targetEdgeBrlen, _targetEdgeInheriProb, _destinationEdge, null, null, _sourceEdge, null, null);
        undo.performOperation();
    }
}
