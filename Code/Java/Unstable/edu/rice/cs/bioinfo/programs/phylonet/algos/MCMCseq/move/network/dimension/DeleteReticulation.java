package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.dimension;

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
public class DeleteReticulation extends DimensionChange {

    private double _logHR;
    private NetNode<NetNodeInfo> _v1, _v2, _v3, _v4, _v5, _v6;
    private double _oldGamma;

    public DeleteReticulation(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _v1 = _v2 = _v3 = _v4 = _v5 = _v6 = null; // reset

        List<NetNode<NetNodeInfo>> netNodes = IterableHelp.toList(_network.getNetwork().getNetworkNodes());
        _v2 = netNodes.get(Randomizer.getRandomInt(netNodes.size()));
        List<NetNode<NetNodeInfo>> parents = IterableHelp.toList(_v2.getParents());
        _v1 = parents.get(Randomizer.getRandomInt(parents.size()));

        if(/*_v1.isRoot() ||*/ _v1.isNetworkNode()) {
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
                    double numRetiEdges = 2 * numRetiNodes;
                    double numEdges = 2 * (_network.getNetwork().getLeafCount() - 1) + 3 * (numRetiNodes - 1);
                    double l1 = (_v3 != null ? _v3.getData().getHeight() : Utils.ROOT_TIME_UPPER_BOUND) - _v4.getData().getHeight();
                    double l2 = _v5.getData().getHeight() - _v6.getData().getHeight();

                    double logPopSize = removeReticulation(_v1, _v2, _v3, _v4, _v5, _v6);

                    _violate = true;
                    // convert to coalescent unit
                    _logHR = Math.log(pad / l1 * _time_scale / l2 * _time_scale
                            * numRetiEdges / numEdges / (numEdges-1.0)) + logPopSize;
                }
            }
        }
        return _logHR;
    }

    @Override
    public void undo() {
        if(_logHR == Utils.INVALID_MOVE) return;
        addReticulation(_v1, _v2, _v3, _v4, _v5, _v6, _oldGamma); // TODO: this is a bug: old parameters may be changed
    }

    @Override
    public String getName() {
        return "Delete-Reticulation";
    }
}
