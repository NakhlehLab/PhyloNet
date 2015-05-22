package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;



import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;

/**
 * Created by yunyu on 10/20/14.
 */
public class ReticulationEdgeDeletion extends NetworkRearrangementOperation{

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
                //_network.addEdge(_sourceEdge);
            }

            //_network.addEdge(_destinationEdge);

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
        //_destinationEdge.Item2.setParentProbability(_destinationEdge.Item1, _destinationEdge.Item2.getParentProbability(_destinationEdge.Item1)/(1-_targetEdgeInheriProb));

        if(_sourceEdge.Item1==null){
            _network.resetRoot(_sourceEdge.Item2);
            _sourceEdge = null;
        }

        return true;
    }



    public void undoOperation(){
        ReticulationEdgeAddition undo = new ReticulationEdgeAddition();
        undo.setParameters(_network, _targetEdge, _targetEdgeBrlen, _targetEdgeInheriProb, _sourceEdge, _sourceEdgeBrlens, _sourceEdgeInheriProbs, _destinationEdge, _destinationEdgeBrlens, _destinationEdgeInheriProbs);
        undo.performOperation();
    }
}
