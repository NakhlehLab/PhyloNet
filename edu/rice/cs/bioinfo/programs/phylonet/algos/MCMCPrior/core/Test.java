package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.core;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.InferTreeWrapper;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityIntegrated;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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

        Network<Double> trueNetwork;
        trueNetwork = Networks.readNetwork("((a:0.1,(b:0.05)#H1:0.05::0.5):0.1,(#H1:0.05::0.5,c:0.1):0.1);");
        trueNetwork = Networks.readNetwork("((a:0.1,b:0.1):0.1,c:0.2);");
        trueNetwork = Networks.readNetwork("(((a:1.0,b:1.0):1.0,c:2.0):1.0,d:3.0);");
        //trueNetwork = Networks.readNetwork("(a:0.1,b:0.1);");
        //trueNetwork = Networks.readNetwork("((((b:1.0,c:1.0)I4:" + y + ")I3#H1:0.1::" + (1 - gamma) + ",a:" + (1.1 + y) + ")I1:" + x + ",(I3#H1:0.1::" + gamma + ",d:" + (1.1 + y) + ")I2:"+ x +")I0;");


        Utils._START_NET = "[0.1]" + trueNetwork.toString();
        System.out.println(Utils._START_NET);

        Utils._CHAIN_LEN = 1500000;
        Utils._BURNIN_LEN = 500000;
        Utils._SAMPLE_FREQUENCY = 500;
        Utils._BL_EXP_PRIOR = true;
        Utils._TIMES_EXP_PRIOR = false;
        Utils.EXP_PARAM = 40.0;
        Utils._CONST_POP_SIZE = false;
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._ESTIMATE_POP_PARAM = false;
        Utils._NETWORK_SIZE_PRIOR = false;
        Utils._SEED = 12345678;
        Utils.DISABLE_TOPOLOGY_MOVES = true;

        //System.out.println(new ExponentialDistribution(1.0).cumulativeProbability(0.1));

        MC3Core run = new MC3Core();
        run.run();

        List<Tree> gts = new ArrayList<>();

        double maxBL = 0.0;
        double avgBL = 0.0;
        double totalB = 0.0;

        for(int i = (int)(Utils._BURNIN_LEN / Utils._SAMPLE_FREQUENCY) ; i < Utils._CHAIN_LEN / Utils._SAMPLE_FREQUENCY ; i++) {
            SimGTInNetworkByMS simgt = new SimGTInNetworkByMS();
            Network<Double> net = Networks.readNetworkWithRootPop(run._networkList.get(i));
            for(NetNode node : net.bfs()) {
                for(Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    maxBL = Math.max(maxBL, node.getParentDistance(parent));
                    totalB += 1;
                    avgBL += node.getParentDistance(parent);
                }
            }

            gts.addAll(simgt.generateGTs(Networks.readNetworkWithRootPop(run._networkList.get(i)), null, 10, "/Users/zhujiafan/Documents/Luay/msdir/ms"));
        }

        avgBL /= totalB;

        System.out.println("MAX BL " + maxBL + " AVG BL" + avgBL);

        List<List<MutableTuple<Tree, Double>>> treeList = new ArrayList<>();
        for(Tree tree : gts) {
            treeList.add(Arrays.asList(new MutableTuple<Tree, Double>(tree, 1.0)));

        }

        InferTreeWrapper inferTreeWrapper = new InferTreeWrapper();
        List<MutableTuple<Tree, Double>> distinctTrees = new ArrayList<>();
        inferTreeWrapper.summarizeData(treeList, null, distinctTrees);

        List<Tree> summarizedGTs = new ArrayList<>();
        for(MutableTuple<Tree, Double> gtTuple : distinctTrees) {
            summarizedGTs.add(gtTuple.Item1);
        }

        GeneTreeProbabilityIntegrated probCalc = new GeneTreeProbabilityIntegrated();
        List<Double> prob = probCalc.calculateGTDistribution(trueNetwork, summarizedGTs, null, 1.0, false);

        double sum = 0.0;
        for(int i = 0 ; i < distinctTrees.size() ; i++) {
            System.out.println(distinctTrees.get(i).Item1.toNewick() + " " + distinctTrees.get(i).Item2 + " " + prob.get(i) * gts.size());
            sum += prob.get(i);
        }

        System.out.println("Sum of prob: " + sum);

    }
}
