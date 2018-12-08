package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.dimension;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.distribution.ExponentialDistribution;

import java.util.List;
import java.util.concurrent.atomic.DoubleAdder;

/**
 * Created by wendingqiao on 3/6/16.
 * Add reticulation
 */
public class AddReticulation extends DimensionChange {

    protected NetNode<NetNodeInfo> _v1, _v2, _v3, _v4, _v5, _v6;

    public AddReticulation(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _v1 = _v2 = _v3 = _v4 = _v5 = _v6 = null; // reset

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = Networks.getAllEdges(_network.getNetwork());
        edges.add(new Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>(_network.getNetwork().getRoot(), null));

        for(Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> tuple : edges) {
            if(!tuple.Item1.isRoot() && !Double.isNaN(tuple.Item1.getRootPopSize())) {
                System.err.println("Not root but has root pop size!");
            }
        }

        int numEdges = edges.size();
        int numRetiNodes = _network.getNetwork().getReticulationCount();
        double rootHeight = _network.getNetwork().getRoot().getData().getHeight();

        Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge1, edge2;
        edge1 = edges.get(Randomizer.getRandomInt(numEdges));
        do {
            edge2 = edges.get(Randomizer.getRandomInt(numEdges));
        } while (edge1 == edge2);

        _v3 = edge1.Item2;
        _v4 = edge1.Item1;
        _v5 = edge2.Item2;
        _v6 = edge2.Item1;

        double rootDelta = -1;
        if(_v3 == null || _v5 == null) rootDelta = nodeHeights.sample();

        double t3 = _v3 != null ? _v3.getData().getHeight() : (rootHeight + rootDelta);
        double t4 = _v4.getData().getHeight();
        double l1 = _v3 != null ? t3 - t4 : _time_scale;
        double t1 = _v3 != null ? t4 + Randomizer.getRandomDouble() * l1 : t4 + rootDelta;

        double t5 = _v5 != null ? _v5.getData().getHeight() : (rootHeight + rootDelta);
        double t6 = _v6.getData().getHeight();
        double l2 = _v5 != null ? t5 - t6 : _time_scale;
        double t2 = _v5 != null ? t6 + Randomizer.getRandomDouble() * l2 : t6 + rootDelta;

        if(l1 <= 0 || l2 <= 0) {
            throw new RuntimeException("Add reticulation error!");
        }

        double gamma = Randomizer.getRandomDouble();
        _v1 = new BniNetNode<>();
        _v1.setData(new NetNodeInfo(t1));
        _v2 = new BniNetNode<>();
        _v2.setData(new NetNodeInfo(t2));

        double logPopSize = (t1 > t2) ? addReticulation(_v1, _v2, _v3, _v4, _v5, _v6, gamma) // v3 v4 t1 v5 v6 t2
                                    : addReticulation(_v2, _v1, _v5, _v6, _v3, _v4, gamma);  // v5 v6 t2 v3 v4 t1

        double pda = numRetiNodes == 0 ? 0.5 : 1.0;
        double numRetiEdges = 2 * (numRetiNodes + 1);

        // convert to the scale of the root population size
        // if only one population size, it equals to convert to coalescent unit
        double hr = pda * l1 / _time_scale * l2 / _time_scale * numEdges * (numEdges-1.0) / numRetiEdges;

        _violate = false;

        double logHasting = Math.log(hr) + logPopSize;
        if(rootDelta > 0) {
            double lambda = 1.0 / nodeHeights.getMean();
            logHasting += -Math.log(lambda) + lambda * rootDelta;
        }

        return logHasting;
    }

    @Override
    public void undo() {
        if(_v1 == null && _v2 == null && _v3 == null && _v4 == null && _v5 == null && _v6 == null) {
            return;
        }

        if(_v1.getData().getHeight() > _v2.getData().getHeight()) {
            removeReticulation(_v1, _v2, _v3, _v4, _v5, _v6);
        } else {
            removeReticulation(_v2, _v1, _v5, _v6, _v3, _v4);
        }
    }

    @Override
    public String getName() {
        return "Add-Reticulation";
    }

}
