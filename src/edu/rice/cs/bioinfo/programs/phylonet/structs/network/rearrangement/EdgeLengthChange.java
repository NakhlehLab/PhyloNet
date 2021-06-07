package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of EdgeParameterChange.
 * It changes the branch length of an edge
 */
public class EdgeLengthChange extends EdgeParameterChange{

    /**
     * This is the main function for changing the branch length of an edge
     *
     * @return  whether the operation is performed successfully
     *          returns false when the resulting network is invalid
     */
    public boolean performOperation(){
        double newBranchLength = _targetEdgeBrlen;
        if(newBranchLength==-1){
            _targetEdgeBrlen = _targetEdge.Item2.getParentDistance(_targetEdge.Item1);
            newBranchLength = drawRandomParameter(_targetEdgeBrlen, 0, Double.POSITIVE_INFINITY);
        }
        _targetEdge.Item2.setParentDistance(_targetEdge.Item1, newBranchLength);
        return true;
    }

    /**
     * This function is to undo the operation
     */
    public void undoOperation(){
        performOperation();
    }
}
