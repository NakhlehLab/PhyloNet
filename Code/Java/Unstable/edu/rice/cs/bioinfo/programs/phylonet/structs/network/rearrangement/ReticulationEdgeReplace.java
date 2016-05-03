package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkRearrangementOperation.
 * It re-places an reticulation edge
 */
public class ReticulationEdgeReplace extends NetworkRearrangementOperation {
    ReticulationEdgeDeletion _deletion;
    ReticulationEdgeAddition _addition;

    /**
     * This is the main function for re-placing an reticulation edge
     *
     * @return  whether the operation is performed successfully
     *          returns false when the resulting network is invalid
     */
    public boolean performOperation(){
        _addition = new ReticulationEdgeAddition();
        _addition.setParameters(_network, null, _sourceEdge, _destinationEdge);
        _addition.performOperation();

        _deletion = new ReticulationEdgeDeletion();
        _deletion.setParameters(_network, _targetEdge, null, null);
        if(_deletion.performOperation()){
            return true;
        }
        else{
            _addition.undoOperation();
            return false;
        }
    }


    /**
     * This function is to undo the operation
     */
    public void undoOperation(){
        _deletion.undoOperation();
        _addition.undoOperation();

    }
}
