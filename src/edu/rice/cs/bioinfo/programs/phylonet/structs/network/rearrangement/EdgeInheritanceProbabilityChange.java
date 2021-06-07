package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of EdgeParameterChange.
 * It changes the inheritance probability of an edge
 */
public class EdgeInheritanceProbabilityChange extends EdgeParameterChange{

    /**
     * This is the main function for changing parameters of an edge
     *
     * @return  whether the operation is performed successfully
     *          returns false when the resulting network is invalid
     */
    public boolean performOperation(){
        double newInheritanceProbability = _targetEdgeInheriProb;
        if(newInheritanceProbability==-1){
            _targetEdgeInheriProb = _targetEdge.Item2.getParentProbability(_targetEdge.Item1);
            newInheritanceProbability = drawRandomParameter(_targetEdgeInheriProb, 0, 1);
        }
        _targetEdge.Item2.setParentProbability(_targetEdge.Item1, newInheritanceProbability);
        NetNode anotherParent = findAnotherParentAndChild(_targetEdge.Item2, _targetEdge.Item1).Item1;
        _targetEdge.Item2.setParentProbability(anotherParent, 1-newInheritanceProbability);
        return true;
    }

    /**
     * This function is to undo the operation
     */
    public void undoOperation(){
        performOperation();
    }
}
