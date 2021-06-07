package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.param;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by wendingqiao on 3/4/16.
 *
 * select a random internal node (not root nor leaf)
 * select a new height between [lower by children & constraints, upper by parents]
 */
public class ChangeTime extends NetworkOperator {

    private NetNode<NetNodeInfo> _target;
    private Double _oldHeight;

    public ChangeTime(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        // reset
        _target = null;
        _oldHeight = null;

        List<NetNode<NetNodeInfo>> nodes = Networks.getInternalNodes(_network.getNetwork());
        do {
            _target = nodes.get(Randomizer.getRandomInt(nodes.size()));
        } while (_target.isRoot());

        double[] bounds = _network.getLowerAndUpperBound(_target);
        // increase height => may violate
        _oldHeight = _target.getData().getHeight();
        if((bounds[0] - _oldHeight) > 0.000001 || (_oldHeight - bounds[1]) > 0.000001) {
            throw new RuntimeException("Invalid height " + _target.getName() + " " + _oldHeight + " "
                    + Arrays.toString(bounds) + "\n" + _network.getNetwork().toString());
        }

        double newHeight = bounds[0] + Randomizer.getRandomDouble() * (bounds[1] - bounds[0]);
        setNodeHeight(_target, newHeight);

        _violate = newHeight > _oldHeight;
        return 0;
    }

    @Override
    public void undo() {
        setNodeHeight(_target, _oldHeight);
    }

    @Override
    public String getName() {
        return "Change-Time";
    }


}
