package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/12/12
 * Time: 4:30 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReticulationEdgeDestinationChange extends NetworkRearrangementOperation{

    public boolean performOperation(){
        if(_sourceEdge == null){
            _sourceEdge = findAnotherParentAndChild(_targetEdge.Item2, _targetEdge.Item1);
            if(_sourceEdge.Item2.hasParent(_sourceEdge.Item1)){
                return false;
            }
        }
        _sourceEdgeBrlens = new double[2];
        _sourceEdgeInheriProbs = new double[2];
        removeNodeFromAnEdge(_targetEdge.Item2, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);
        if(_destinationEdgeBrlens==null){
            _destinationEdgeBrlens = new double[2];
            _destinationEdgeInheriProbs = new double[2];
            randomlyPartitionAnEdge(_destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
        }
        addNodeToAnEdge(_targetEdge.Item2, _destinationEdge,_destinationEdgeBrlens, _destinationEdgeInheriProbs);
        _targetEdge.Item2.setParentProbability(_destinationEdge.Item1, _targetEdge.Item2.getParentProbability(_targetEdge.Item1)==NetNode.NO_PROBABILITY? NetNode.NO_PROBABILITY: 1-_targetEdge.Item2.getParentProbability(_targetEdge.Item1));
        return true;
    }

    public void undoOperation(){
        ReticulationEdgeDestinationChange undo = new ReticulationEdgeDestinationChange();
        undo.setParameters(_network, _targetEdge, -1, -1, _destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);
        undo.performOperation();
    }
    



}
