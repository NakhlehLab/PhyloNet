package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/7/12
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReticulationEdgeReplace extends NetworkRearrangementOperation {
    ReticulationEdgeDeletion _deletion;
    ReticulationEdgeAddition _addition;

    /*
    public boolean performOperation(){
        _deletion = new ReticulationEdgeDeletion();
        _deletion.setParameters(_network, _targetEdge, null, null);
        if(_deletion.performOperation()){
            _addition = new ReticulationEdgeAddition();
            _addition.setParameters(_network, null, _sourceEdge, _destinationEdge);
            _addition.performOperation();
            return true;
        }
        else{
            return false;
        }
    }
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


    public void undoOperation(){
        _deletion.undoOperation();
        _addition.undoOperation();

    }
}
