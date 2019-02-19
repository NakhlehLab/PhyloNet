package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.dimension;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.List;

/**
 * Created by wendingqiao on 3/6/16.
 * Delete reticulation
 */
public class DeleteReticulationFromBackbone extends DeleteReticulation {

    private Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> retiEdge = null;

    public DeleteReticulationFromBackbone(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _v1 = _v2 = _v3 = _v4 = _v5 = _v6 = null; // reset

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> retiEdges = _network.getRetiEdges();
        retiEdge = retiEdges.get(Randomizer.getRandomInt(retiEdges.size()));
        _v2 = retiEdge.Item1;
        _v1 = retiEdge.Item2;

        if(/*_v1.isRoot() ||*/ _v1.isNetworkNode() || !_v2.hasParent(_v1)) {
            _violate = false;
            _logHR = Utils.INVALID_MOVE;
        } else {
            _v3 = _v1.isRoot() ? null : _v1.getParents().iterator().next();
            _v4 = Networks.getOtherChild(_v1, _v2);
            if(_v4.hasParent(_v3)) {
                _violate = false;
                _logHR = Utils.INVALID_MOVE;
            } else {
                _v5 = Networks.getOtherParent(_v2, _v1);
                _v6 = _v2.getChildren().iterator().next();
                if(_v6.hasParent(_v5) || (_v3 == _v5 && _v4 == _v6)) {
                    _violate = false;
                    _logHR = Utils.INVALID_MOVE;
                } else {
                    _oldGamma = _v2.getParentProbability(_v1);
                    int numRetiNodes = _network.getNetwork().getReticulationCount();
                    double pad = numRetiNodes == 1 ? 2.0 : 1.0;
                    double numEdges = 2 * (_network.getNetwork().getLeafCount() - 1) + 3 * (numRetiNodes - 1);
                    double l1 = _v3 != null ? (_v3.getData().getHeight() - _v4.getData().getHeight()) : _time_scale ;
                    double l2 = _v5.getData().getHeight() - _v6.getData().getHeight();

                    double logPopSize = removeReticulation(_v1, _v2, _v3, _v4, _v5, _v6);
                    _network.getRetiEdges().remove(retiEdge);

                    _violate = true;
                    // convert to coalescent unit
                    _logHR = Math.log(pad / l1 * _time_scale / l2 * _time_scale
                            * numRetiNodes / numEdges / (numEdges-1.0)) + logPopSize;

                    if(_v3 == null) {
                        double lambda = 1.0 / nodeHeights.getMean();
                        double rootDelta = _v1.getData().getHeight() - _v4.getData().getHeight();
                        _logHR += Math.log(lambda) - lambda * rootDelta;
                    }
                }
            }
        }
        if(_logHR == Utils.INVALID_MOVE) {
            retiEdge = null;
        }

        return _logHR;
    }

    @Override
    public void undo() {
        super.undo();
        if(retiEdge != null) {
            _network.getRetiEdges().add(retiEdge);
            retiEdge = null;
        }
    }

    @Override
    public String getName() {
        return "Delete-Reticulation-From-Backbone";
    }
}
