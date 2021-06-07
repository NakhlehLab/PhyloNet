package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal;

/**
 * Change branch length of current network
 * Created by dw20 on 5/31/15.
 * parameters: _targetEdge
 */
public class ChangeLength extends NetworkGTTOperation {

    private double hastingsRatio = Double.MIN_VALUE;
    private double prevLength = -1;
    private double windowSize = 0.10; // [-ws/2, +ws/2]

    public double getLogHR() {
        return hastingsRatio;
    }

    public boolean performOperation() {
        if(_targetEdge == null) return false;
        prevLength = _targetEdge.Item2.getParentDistance(_targetEdge.Item1); // get prev length
//        hastingsRatio = 0.0;
//        double newBranch = prevLength + windowSize * (_random.nextDouble() - 0.5) * 2.0;
//        if(newBranch < 0.0) {
//            newBranch = prevLength + windowSize * (_random.nextDouble() - 0.5) * 2.0;
//        }
        hastingsRatio = windowSize * (_random.nextDouble() - 0.5); // set hasting ratio to e^(u-0.5)  -> log HR = u-0.5
        double newBranch = prevLength * Math.exp(hastingsRatio);
        _targetEdge.Item2.setParentDistance(_targetEdge.Item1, newBranch); // set new length
        return true;
    }

    public void undoOperation() {
        _targetEdge.Item2.setParentDistance(_targetEdge.Item1, prevLength);
        hastingsRatio = Double.MIN_VALUE;
    }
}
