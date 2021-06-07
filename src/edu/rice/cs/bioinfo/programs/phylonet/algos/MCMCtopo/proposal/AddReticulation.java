package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal;


import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

/**
 * Add reticulation to current network.
 *
 * Created by dw20 on 5/31/15.
 * parameters: _sourceEdge _destEdge
 */
public class AddReticulation extends NetworkGTTOperation {

    private double hastingsRatio = Double.MIN_VALUE;

    public boolean performOperation(){

        if(_sourceEdge == null) return false;

        int reti = _network.getReticulationCount();
        double re = 2.0 * reti + 2.0; // number of reti edges after addition
        double ne = 2.0 * _network.getLeafCount() + 3.0 * reti - 2.0; // number of edges in current network

        NetNode newEdgeTail = new BniNetNode();
        NetNode newEdgeHead = new BniNetNode();

        if(_targetEdge != null) {
            newEdgeTail = _targetEdge.Item1;
            newEdgeHead = _targetEdge.Item2;
        }
        double srcEdgeLen = 0, destEdgeLen = 0;

        if(_sourceEdge!=null){
            if(_sourceEdgeBrlens==null) {
                _sourceEdgeBrlens = new double[2];
                _sourceEdgeInheriProbs = new double[2];
                randomlyPartitionAnEdge(_sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);
            }
            srcEdgeLen = _sourceEdge.Item2.getParentDistance(_sourceEdge.Item1);

            addNodeToAnEdge(newEdgeTail, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);
        }

        if(_destinationEdgeBrlens==null) {
            _destinationEdgeBrlens = new double[2];
            _destinationEdgeInheriProbs = new double[2];
            randomlyPartitionAnEdge(_destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
        }
        destEdgeLen = _destinationEdge.Item2.getParentDistance(_destinationEdge.Item1);
        addNodeToAnEdge(newEdgeHead, _destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);


        _targetEdge = new Tuple<>(newEdgeTail, newEdgeHead);
        if(_targetEdgeBrlen==-1) {
            double w1 = _random.nextDouble(); // w1 = 1.0 - w1
            _targetEdgeBrlen = -Math.log(w1); // propose a new branch length
            _targetEdgeInheriProb = _random.nextDouble(); // propose a new probability from U(0,1)
            hastingsRatio = ne * (ne-1.0) * srcEdgeLen * destEdgeLen / (re * Math.exp(-_targetEdgeBrlen));
//            System.out.println("Hastings-Ratio-Add-Reticulation: " + srcEdgeLen + " " + destEdgeLen + " " + _targetEdgeBrlen);
        }
        newEdgeTail.adoptChild(newEdgeHead, _targetEdgeBrlen);
        newEdgeHead.setParentProbability(newEdgeTail, _targetEdgeInheriProb);
        newEdgeHead.setParentProbability(_destinationEdge.Item1, 1-_targetEdgeInheriProb);

        Networks.autoLabelNodes(_network);

        return true;
    }

    public void undoOperation(){
        DeleteReticulation undo = new DeleteReticulation();
        undo.setParameters(_network, _targetEdge, -1, -1,
                null, _sourceEdgeBrlens, _sourceEdgeInheriProbs,
                null, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
        undo.performOperation();
        hastingsRatio = Double.MIN_VALUE;
    }

    public double getLogHR() {
        return hastingsRatio == Double.MIN_VALUE ? Double.MIN_VALUE : Math.log(hastingsRatio);
    }

}
