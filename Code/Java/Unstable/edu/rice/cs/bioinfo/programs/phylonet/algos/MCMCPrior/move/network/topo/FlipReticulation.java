package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.topo;

import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.distribution.SpeciesNetPriorDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 3/7/16.
 *
 * select a random network node
 * select a random reticulation edge [node, parent]
 * if not removable (parent is also a network node) => nullify the move
 * flip the edge
 */
public class FlipReticulation extends NetworkOperator {

    private double _logHR;
    private NetNode<NetNodeInfo> _v1, _v2; // original child/parent
    private double _t1, _t2;

    public FlipReticulation(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _v1 = _v2 = null; // reset

        List<NetNode<NetNodeInfo>> netNodes = IterableHelp.toList(_network.getNetwork().getNetworkNodes());
        _v2 = netNodes.get(Randomizer.getRandomInt(netNodes.size()));
        List<NetNode<NetNodeInfo>> parents = IterableHelp.toList(_v2.getParents());
        _v1 = parents.get(Randomizer.getRandomInt(parents.size()));

        if(_v1.isNetworkNode()) {
            _logHR = Utils.INVALID_MOVE;
            _violate = false;
        } else {

            NetNode<NetNodeInfo> v4, v5;
            v4 = Networks.getOtherChild(_v1, _v2);
            v5 = Networks.getOtherParent(_v2, _v1);
            double t4 = v4.getData().getHeight();
            double t5 = v5.getData().getHeight();

            if(t4 >= t5) {
                _logHR = Utils.INVALID_MOVE;
                _violate = false;
            } else {
                _t1 = _v1.getData().getHeight();
                _t2 = _v2.getData().getHeight();
                // v1 = parent, v2 = target
                NetNode<NetNodeInfo> v3, v6;
                v3 = _v1.getParents().iterator().next();
                v6 = _v2.getChildren().iterator().next();
                double t3 = v3.getData().getHeight();
                double t6 = v6.getData().getHeight();

                double th = Math.min(t3, t5);
                double tl = Math.max(t6, t4);
                double t1Prime = t4 + Randomizer.getRandomDouble() * (th - t4);
                double t2Prime = tl + Randomizer.getRandomDouble() * (t5 - tl);
                if(t1Prime > t2Prime) {
                    double temp = t1Prime;
                    t1Prime = t2Prime;
                    t2Prime = temp;
                }
                flip(_v1, _v2, t1Prime, t2Prime);
                double logPrimeP = - Math.log(t5 - tl) - Math.log(th - t4);
                if(tl < t1Prime && t1Prime < th && tl < t2Prime && t2Prime < th) {
                    logPrimeP += Math.log(2.0);
                }
                double logPPrime = - Math.log(t3 - tl) - Math.log(th - t6);
                if(tl < _t1 && _t1 < th && tl < _t2 && _t2 < th) {
                    logPPrime += Math.log(2.0);
                }

                _logHR = logPPrime - logPrimeP;
                _violate = true;
            }
        }
        return _logHR;
    }

    @Override
    public void undo() {
        if(_logHR == Utils.INVALID_MOVE) return;
        flip(_v2, _v1, _t2, _t1);
    }

    @Override
    public String getName() {
        return "Flip-Reticulation";
    }

    // child is prev parent, current child
    // parent is prev child, current parent
    // cHeight is new child height
    private void flip(NetNode<NetNodeInfo> child, NetNode<NetNodeInfo> parent, double cHeight, double pHeight) {
        setNodeHeight(child, cHeight);
        setNodeHeight(parent, pHeight);

        double[] params = getParameters(child, parent);

        child.removeChild(parent);

        adopt(parent, child, params);

        child.setParentProbability(Networks.getOtherParent(child, parent), 1.0 - params[0]);

        parent.setParentProbability(parent.getParents().iterator().next(), Double.NaN);
    }



}
