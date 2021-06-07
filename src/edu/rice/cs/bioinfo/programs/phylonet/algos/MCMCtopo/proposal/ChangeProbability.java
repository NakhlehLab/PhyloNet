package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.proposal;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.EdgeParameterChange;

/**
 * Change a inheritance probability of current network
 * Created by dw20 on 06/07/15.
 * parameters: _targetEdge
 */
public class ChangeProbability extends NetworkGTTOperation{

    private double windowSize = 0.1; //[-ws/2, +ws/2]
    private double prevProb = -1;

    public double getLogHR() {
        return 0.0;  // always return log(1.0) = 0.0;
    }

    public boolean performOperation(){
        if(_targetEdge == null) return false;
        prevProb = _targetEdge.Item2.getParentProbability(_targetEdge.Item1);
        // set new prob
        double newProb = prevProb + (_random.nextDouble()-0.50) * windowSize;
        if(newProb < 0.0) {
            newProb = -newProb;
        } else if(newProb > 1.0) {
            newProb = 2.0 - newProb;
        }
        setProb(newProb);
        return true;
    }

    public void undoOperation(){
        setProb(prevProb);
    }

    private void setProb(double prob) {
        _targetEdge.Item2.setParentProbability(_targetEdge.Item1, prob);
        NetNode anotherParent = findAnotherParentAndChild(_targetEdge.Item2, _targetEdge.Item1).Item1;
        _targetEdge.Item2.setParentProbability(anotherParent, 1.0 - prob);
    }
}
