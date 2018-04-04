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
 * Scales the time of the root by a random factor.
 */
public class ScaleRootTime extends NetworkOperator {

    private double _scaleFactor = 0.92;
    private double _upperLimit = 1.0 - 1e-6;
    private double _lowerLimit = 0.85;

    public double scale;
    private double _logHR;
    private Double _oldHeight;

    public ScaleRootTime(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        // reset
        _oldHeight = null;

        scale = getScaler(); // scale > 1 => increase height => may violate

        NetNode<NetNodeInfo> root = _network.getNetwork().getRoot();
        _oldHeight = root.getData().getHeight();
        double newHeight = _oldHeight * scale;
        double[] bounds = _network.getLowerAndUpperBound(root);

        if(newHeight < bounds[0] || newHeight > Utils.NET_MAX_HEIGHT) {
            _logHR = Utils.INVALID_MOVE;
            _violate = false;
        } else {
            setNodeHeight(root, newHeight);
            _logHR = -Math.log(scale);
            _violate = scale > 1.0;
        }
        return _logHR;
    }

    @Override
    public void undo() {
        if (_logHR == Utils.INVALID_MOVE) return;
        setNodeHeight(_network.getNetwork().getRoot(), _oldHeight);
    }

    @Override
    public String getName() {
        return "Scale-Root-Time";
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

}