package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal;

/**
 * Delete a reticulation edge from current network
 * Created by dw20 on 5/31/15.
 * parameters: _targetEdge
 */
public class DeleteReticulation extends NetworkGTTOperation {

    private double hastingsRatio = Double.MIN_VALUE;
    public double srcDistLen = -1.0, destDistLen = -1.0;

    public boolean performOperation(){
        int reti = _network.getReticulationCount();
        double re = 2.0 * reti; // number of reti edges
        double ne = 2.0 * _network.getLeafCount() + 3.0 * reti - 5.0;  // number of edges of the proposed network (after deletion)

        if(_sourceEdge == null){
            _sourceEdge = findParentAndAnotherChild(_targetEdge.Item1, _targetEdge.Item2);
            _destinationEdge = findAnotherParentAndChild(_targetEdge.Item2, _targetEdge.Item1);
            if(_destinationEdge.Item2.hasParent(_destinationEdge.Item1)){
                return false;
            }
            if(_sourceEdge.Item1!=null){
                if(_sourceEdge.Item2.hasParent(_sourceEdge.Item1) || _destinationEdge.equals(_sourceEdge)){
                    return false;
                }
            }
            if(_sourceEdge.Item1==null){
                return false;
            }
        }


        _targetEdgeBrlen = _targetEdge.Item2.getParentDistance(_targetEdge.Item1);
        _targetEdgeInheriProb = _targetEdge.Item2.getParentProbability(_targetEdge.Item1);
        _targetEdge.Item2.removeChild(_targetEdge.Item1);

        _destinationEdgeBrlens = new double[2];
        _destinationEdgeInheriProbs = new double[2];
        removeNodeFromAnEdge(_targetEdge.Item2, _destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);

        _sourceEdgeBrlens = new double[2];
        _sourceEdgeInheriProbs = new double[2];
        removeNodeFromAnEdge(_targetEdge.Item1, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);
        //_destinationEdge.Item2.setParentProbability(_destinationEdge.Item1, _destinationEdge.Item2.getParentProbability(_destinationEdge.Item1)/(1-_targetEdgeInheriProb));

        srcDistLen  = _sourceEdge.Item2.getParentDistance(_sourceEdge.Item1);
        destDistLen = _destinationEdge.Item2.getParentDistance(_destinationEdge.Item1);

        hastingsRatio = re * Math.exp(-_targetEdgeBrlen) / ( ne * (ne-1.0) * srcDistLen * destDistLen );
//        System.out.println("Hastings-Ratio-Delete-Reticulation: " + srcDistLen + " " + destDistLen + " " + _targetEdgeBrlen);

        return true;
    }



    public void undoOperation(){
        AddReticulation undo = new AddReticulation();
        undo.setParameters(_network, _targetEdge, _targetEdgeBrlen, _targetEdgeInheriProb,
                _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs,
                _destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
        undo.performOperation();
        hastingsRatio = Double.MIN_VALUE;
    }

    public double getLogHR() {
        return hastingsRatio == Double.MIN_VALUE ? Double.MIN_VALUE : Math.log(hastingsRatio);
    }
}
