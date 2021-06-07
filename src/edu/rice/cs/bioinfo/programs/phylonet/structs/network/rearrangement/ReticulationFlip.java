package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkRearrangementOperation.
 * It flips the direction of an reticulation edge
 */
public class ReticulationFlip extends NetworkRearrangementOperation {


    /**
     * This is the main function for flipping the direction an reticulation edge
     *
     * @return  whether the operation is performed successfully
     *          returns false when the resulting network is invalid
     */
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


    /**
     * This function is to undo the operation
     */
    public void undoOperation(){
        ReticulationFlip undo = new ReticulationFlip();
        undo.setParameters(_network, _targetEdge, _targetEdgeBrlen, _targetEdgeInheriProb, _destinationEdge, null, null, _sourceEdge, null, null);
        undo.performOperation();
    }

}
