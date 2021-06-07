package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * Move the tail of a random edge from current network
 * Created by dw20 on 5/31/15.
 * parameters: _targetEdge _destinationEdge
 */
public class MoveTail extends NetworkGTTOperation {

    private double hastingRatio = Double.MIN_VALUE;
    private double rootEdgeLength = 1.0;

    public boolean performOperation(){
        if(_sourceEdge == null){
            _sourceEdge = findParentAndAnotherChild(_targetEdge.Item1, _targetEdge.Item2);
            _sourceEdgeBrlens = new double[2];
            _sourceEdgeInheriProbs = new double[2];
        }

        if (_sourceEdge.Item1 != null && _sourceEdge.Item2.hasParent(_sourceEdge.Item1)) {
            return false;
        }
        removeNodeFromAnEdge(_targetEdge.Item1, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);

        if(_sourceEdge.Item1 == null){
            if(_destinationEdge == null) return false;

            rootEdgeLength = _sourceEdge.Item2.getParentDistance(_targetEdge.Item1);

            _network.resetRoot(_sourceEdge.Item2);
            _sourceEdge = null;
        }

        if(_destinationEdge !=null ){
            if(_sourceEdge != null) {
                // dest edge / source edge
                hastingRatio = _destinationEdge.Item2.getParentDistance(_destinationEdge.Item1) /
                        _sourceEdge.Item2.getParentDistance(_sourceEdge.Item1);
            } else {
                // dest edge * deleted root edge
                hastingRatio = _destinationEdge.Item2.getParentDistance(_destinationEdge.Item1) *
                        1.0 * Math.exp( -1.0 * rootEdgeLength );
            }

            if(_destinationEdgeBrlens==null){
                _destinationEdgeBrlens = new double[2];
                _destinationEdgeInheriProbs = new double[2];
                randomlyPartitionAnEdge(_destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
            }
            addNodeToAnEdge(_targetEdge.Item1, _destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);

        }
        else{
            rootEdgeLength = -1.0 * Math.log( 1.0 - _random.nextDouble() );
            hastingRatio = 1.0 / _sourceEdge.Item2.getParentDistance(_sourceEdge.Item1) /
                    Math.exp( -1.0 * rootEdgeLength );

            if(_destinationEdgeBrlens == null) _destinationEdgeBrlens = _sourceEdgeBrlens;
            _destinationEdgeBrlens[1] = rootEdgeLength;

            if(_destinationEdgeInheriProbs == null) {
                _destinationEdgeInheriProbs = new double[]{ NetNode.NO_PROBABILITY, NetNode.NO_PROBABILITY };
            }
            _targetEdge.Item1.adoptChild(_network.getRoot(), _destinationEdgeBrlens[1]);
            _network.getRoot().setParentProbability(_targetEdge.Item1, _destinationEdgeInheriProbs[1]);
            _network.resetRoot(_targetEdge.Item1);
        }
        return true;

    }

    public void undoOperation(){
        setParameters(_network, _targetEdge, -1, -1, null, null, null, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);
        performOperation();
    }

    public double getLogHR() {
        return hastingRatio == Double.MIN_VALUE ? Double.MIN_VALUE : Math.log( hastingRatio );
    }
}
