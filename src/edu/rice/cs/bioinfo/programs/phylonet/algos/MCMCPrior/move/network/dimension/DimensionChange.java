package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.dimension;

import edu.rice.cs.bioinfo.library.programming.Tuple;
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
 * Created by wendingqiao on 3/22/16.
 * Abstract operator that changes dimension of network
 */
public abstract class DimensionChange extends NetworkOperator {

    public DimensionChange(UltrametricNetwork net) {
        super(net);
    }

    protected double removeReticulation(NetNode<NetNodeInfo> v1, NetNode<NetNodeInfo> v2,
                                        NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                                        NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6) {

        double[] paramV1V4 = getParameters(v1, v4);
        double[] paramV2V6 = getParameters(v2, v6);
        double a = Math.min(paramV1V4[1], paramV2V6[1]);
        double b = Math.max(paramV1V4[1], paramV2V6[1]);
        double logProb = Utils.varyPopSizeAcrossBranches() ?
                getLogPopSizeProb(v1.getParentSupport(v3), a, b) +
                getLogPopSizeProb(v2.getParentSupport(v1), a, b) +
                getLogPopSizeProb(v2.getParentSupport(v5), a, b) : 0;

        v3.removeChild(v1);
        v1.removeChild(v4);
        v5.removeChild(v2);
        v2.removeChild(v6);
        v1.removeChild(v2);

        adopt(v3, v4, paramV1V4);
        adopt(v5, v6, paramV2V6);
        return logProb;
    }

    protected double addReticulation(NetNode<NetNodeInfo> v1, NetNode<NetNodeInfo> v2,
                                     NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                                     NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6,
                                     double gamma) {

        double[] paramV3V4 = getParameters(v3, v4);
        double[] paramV5V6 = getParameters(v5, v6);

        double a = Math.min(paramV3V4[1], paramV5V6[1]);
        double b = Math.max(paramV3V4[1], paramV5V6[1]);

        Tuple<Double, Double> popSizeV3V1 = sample(a, b);
        Tuple<Double, Double> popSizeV5V2 = sample(a, b);
        Tuple<Double, Double> popSizeV1V2 = sample(a, b);

        v3.removeChild(v4);
        v5.removeChild(v6);

        adopt(v3, v1, new double[] {NetNode.NO_PROBABILITY, popSizeV3V1.Item1} );
        adopt(v1, v4, paramV3V4);
        adopt(v5, v2, new double[] {1.0-gamma, popSizeV5V2.Item1} );
        adopt(v2, v6, paramV5V6);

        adopt(v1, v2, new double[] {gamma, popSizeV1V2.Item1} );
        return Math.log(popSizeV3V1.Item2) + Math.log(popSizeV5V2.Item2) + Math.log(popSizeV1V2.Item2);
    }

    protected void undoDeleteReticulation(NetNode<NetNodeInfo> v1, NetNode<NetNodeInfo> v2,
                                          NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                                          NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6,
                                          double gamma, double popSizeV3V1, double popSizeV5V2, double popSizeV1V2) {

        double[] paramV3V4 = getParameters(v3, v4);
        double[] paramV5V6 = getParameters(v5, v6);

        double a = Math.min(paramV3V4[1], paramV5V6[1]);
        double b = Math.max(paramV3V4[1], paramV5V6[1]);

        v3.removeChild(v4);
        v5.removeChild(v6);

        adopt(v3, v1, new double[] {NetNode.NO_PROBABILITY, popSizeV3V1} );
        adopt(v1, v4, paramV3V4);
        adopt(v5, v2, new double[] {1.0-gamma, popSizeV5V2} );
        adopt(v2, v6, paramV5V6);

        adopt(v1, v2, new double[] {gamma, popSizeV1V2} );
    }

    // generate sample & partial derivative from sampling distribution
    private Tuple<Double, Double> sample(double a, double b) {
        boolean setPopSize = Utils.varyPopSizeAcrossBranches();
        if(!setPopSize) {
            return new Tuple<>(Double.NEGATIVE_INFINITY, 1.0);
        }
        double c = 1.5 / (2 * b - a);
        double u = Randomizer.getRandomDouble();
        if(u < c * a / 2.0) {
            return new Tuple<>(Math.sqrt(2.0 * a * u / c), Math.sqrt(a / 2.0 / c / u));
        } else if (u <= 0.75) {
            return new Tuple<>(u / c + 0.5 * a, 1.0 / c);
        } else {
            return new Tuple<>(-0.25 / c * Math.log(4.0 * (1.0 - u)) + b, 0.25 / c / (1.0 - u));
        }
    }

    private double getLogPopSizeProb(double x, double a, double b) {
        double c = 1.5 / (2 * b - a);
        if(x < a) {
            return Math.log(c / a * x);
        } else if(x <= b) {
            return Math.log(c);
        } else {
            return Math.log(c) - 4.0 * c * (x - b);
        }
    }



}
