package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;


import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;

import java.util.Arrays;

/**
 * Created by yunyu on 10/20/14.
 */
public class EdgeSourceChange extends NetworkRearrangementOperation{

    public boolean performOperation(){
        if(_sourceEdge == null){
            _sourceEdge = findParentAndAnotherChild(_targetEdge.Item1, _targetEdge.Item2);
            _sourceEdgeBrlens = new double[2];
            _sourceEdgeInheriProbs = new double[2];
        }

        if (_sourceEdge.Item1 != null && _sourceEdge.Item2.hasParent(_sourceEdge.Item1)) {
            return false;
        }
        removeNodeFromAnEdge(_targetEdge.Item1, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);

        if(_sourceEdge.Item1 == null){
            _network.resetRoot(_sourceEdge.Item2);
            _sourceEdge = null;
        }

        if(_destinationEdge !=null ){
            if(_destinationEdgeBrlens==null){
                _destinationEdgeBrlens = new double[2];
                _destinationEdgeInheriProbs = new double[2];
                randomlyPartitionAnEdge(_destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
            }
            addNodeToAnEdge(_targetEdge.Item1, _destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
        }
        else{
            _targetEdge.Item1.adoptChild(_network.getRoot(), _destinationEdgeBrlens[1]);
            _network.getRoot().setParentProbability(_targetEdge.Item1, _destinationEdgeInheriProbs[1]);
            _network.resetRoot(_targetEdge.Item1);
        }
        return true;

    }

    public void undoOperation(){
        setParameters(_network, _targetEdge, -1, -1, null, null, null, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);
        performOperation();
    }
}
