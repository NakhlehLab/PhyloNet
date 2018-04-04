package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.core;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 4/2/18
 * Time: 10:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class Test {
    public static void main(String[] args) {
        generateNetworks();

    }

    public static void generateNetworks() {
        double gamma = 1 - 0.3;

        double y = 1.5;
        double x = 0.5;

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork("((((b:1.0,c:1.0)I4:" + y + ")I3#H1:0.1::" + (1 - gamma) + ",a:" + (1.1 + y) + ")I1:" + x + ",(I3#H1:0.1::" + gamma + ",d:" + (1.1 + y) + ")I2:"+ x +")I0;");

        Utils._START_NET = "[0.1]" + trueNetwork.toString();
        System.out.println(Utils._START_NET);

        Utils._CHAIN_LEN = 1500000;
        Utils._BURNIN_LEN = 500000;
        Utils._SAMPLE_FREQUENCY = 500;
        Utils._BL_EXP_PRIOR = true;
        Utils._TIMES_EXP_PRIOR = false;
        Utils.EXP_PARAM = 1.0;
        Utils._CONST_POP_SIZE = false;
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._ESTIMATE_POP_PARAM = false;
        Utils._NETWORK_SIZE_PRIOR = false;
        Utils._SEED = 12345678;
        Utils._TIME_WINDOW_SIZE = 0.04;
        Utils._POP_SIZE_WINDOW_SIZE = 0.04;
        Utils.DISABLE_TOPOLOGY_MOVES = true;

        //System.out.println(new ExponentialDistribution(1.0).cumulativeProbability(0.1));

        MC3Core run = new MC3Core();
        run.run();

    }
}
