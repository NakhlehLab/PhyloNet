package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

/**
 * Created by yunyu on 10/29/14.
 */
public class EdgeLengthChange extends EdgeParameterChange{
    public boolean performOperation(){
        double newBranchLength = _targetEdgeBrlen;
        if(newBranchLength==-1){
            _targetEdgeBrlen = _targetEdge.Item2.getParentDistance(_targetEdge.Item1);
            newBranchLength = drawRandomParameter(_targetEdgeBrlen, 0, Double.POSITIVE_INFINITY);
        }
        _targetEdge.Item2.setParentDistance(_targetEdge.Item1, newBranchLength);
        return true;
    }

    public void undoOperation(){
        performOperation();
    }
}
