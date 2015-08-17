package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;


import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

/**
 * Created by yunyu on 10/20/14.
 */
public class ReticulationEdgeAddition extends NetworkRearrangementOperation{

    public boolean performOperation(){

        if(_targetEdge == null) {
            _targetEdge = new Tuple<>((NetNode)new BniNetNode(), (NetNode)new BniNetNode());
        }
        NetNode newEdgeTail = _targetEdge.Item1;

        if(_sourceEdge!=null){
            if(_sourceEdgeBrlens==null) {
                _sourceEdgeBrlens = new double[2];
                _sourceEdgeInheriProbs = new double[2];
                randomlyPartitionAnEdge(_sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);
            }
            addNodeToAnEdge(newEdgeTail, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs);

        }else{
            newEdgeTail.adoptChild(_network.getRoot(), _sourceEdgeBrlens[1]);
            _network.getRoot().setParentProbability(newEdgeTail, _sourceEdgeInheriProbs[1]);
            _network.resetRoot(newEdgeTail);
        }

        NetNode newEdgeHead = _targetEdge.Item2;
        if(_destinationEdgeBrlens==null) {
            _destinationEdgeBrlens = new double[2];
            _destinationEdgeInheriProbs = new double[2];
            randomlyPartitionAnEdge(_destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
        }
        addNodeToAnEdge(newEdgeHead, _destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);

        if(_targetEdgeBrlen==-1) {
            _targetEdgeBrlen = 1;
            _targetEdgeInheriProb = 0.5;
        }
        newEdgeTail.adoptChild(newEdgeHead, _targetEdgeBrlen);
        newEdgeHead.setParentProbability(newEdgeTail, _targetEdgeInheriProb);
        newEdgeHead.setParentProbability(_destinationEdge.Item1,_targetEdgeInheriProb==NetNode.NO_PROBABILITY? NetNode.NO_PROBABILITY: 1-_targetEdgeInheriProb);
        Networks.autoLabelNodes(_network);
        return true;
    }

    public void undoOperation(){
        ReticulationEdgeDeletion undo = new ReticulationEdgeDeletion();
        undo.setParameters(_network, _targetEdge, -1, -1, null, _sourceEdgeBrlens, _sourceEdgeInheriProbs, null, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
        undo.performOperation();
    }
}
