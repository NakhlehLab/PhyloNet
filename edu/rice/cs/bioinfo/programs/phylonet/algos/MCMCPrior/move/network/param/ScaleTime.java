package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.param;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 3/4/16.
 * Scale the time of a random selected node by a random factor
 */
public class ScaleTime extends NetworkOperator {

    private double _scaleFactor = 0.95;
    private double _upperLimit = 1.0 - 1e-6;
    private double _lowerLimit = 0.9;

    public Double scale;

    public ScaleTime(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        scale = null; // reset
        scale = getScaler(); // scale > 1 => increase height => may violate
        scaleHeight(scale, false);
        _violate = scale > 1.0;
        return Math.log(scale) * (_network.getInternalNodeCount() - 2);
    }

    @Override
    public void undo() {
        scaleHeight(1.0 / scale, true);
    }

    @Override
    public String getName() {
        return "Scale-Time";
    }

    @Override
    public void optimize(double logAlpha) {
        double delta = Utils.calcDelta(this, logAlpha);
        delta += Math.log(1.0 / _scaleFactor - 1.0);
        setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
    }

    private void setCoercableParameterValue(final double value) {
        _scaleFactor = Math.max(Math.min(value, _upperLimit), _lowerLimit);
    }

    private double getScaler() {
        return _scaleFactor + Randomizer.getRandomDouble() * (1.0 / _scaleFactor - _scaleFactor);
    }

    private void scaleHeight(double scale, boolean undo) {
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_network.getNetwork())) {
            if(node.isLeaf()) continue;
            if(undo) {
                node.getData().recoverHeight();
            } else {
                node.getData().storeHeight();
                node.getData().setHeight(node.getData().getHeight() * scale);
            }
            for(NetNode<NetNodeInfo> child : node.getChildren()) {
                child.setParentDistance(node, node.getData().getHeight() - child.getData().getHeight());
            }
        }
    }


}
