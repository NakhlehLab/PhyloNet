package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;



import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkRearrangementOperation.
 * It removes an reticulation edge from the network
 */
public class ReticulationEdgeDeletion extends NetworkRearrangementOperation{

    /**
     * This is the main function for deleting an reticulation edge
     *
     * @return  whether the operation is performed successfully
     *          returns false when the resulting network is invalid
     */
    public boolean performOperation(){
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

        if(_sourceEdge.Item1==null){
            _network.resetRoot(_sourceEdge.Item2);
            _sourceEdge = null;
        }

        return true;
    }


    /**
     * This function is to undo the operation
     */
    public void undoOperation(){
        ReticulationEdgeAddition undo = new ReticulationEdgeAddition();
        undo.setParameters(_network, _targetEdge, _targetEdgeBrlen, _targetEdgeInheriProb, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs, _destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
        undo.performOperation();
    }
}
