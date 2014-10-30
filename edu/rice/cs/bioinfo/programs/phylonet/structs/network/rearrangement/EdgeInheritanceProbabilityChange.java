package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * Created by yunyu on 10/29/14.
 */
public class EdgeInheritanceProbabilityChange extends EdgeParameterChange{

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

    public void undoOperation(){
        performOperation();
    }
}
