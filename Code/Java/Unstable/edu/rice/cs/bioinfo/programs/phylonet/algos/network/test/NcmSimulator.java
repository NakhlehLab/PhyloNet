package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.core.MC3Core;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.*;
import java.util.concurrent.Exchanger;

/**
 * Created by hunter on 7/12/18.
 */
public class NcmSimulator {
    public double _prior_bl_mean = 1.0;
    public double _prior_bl_max = Double.POSITIVE_INFINITY;

    public static void main(String[] args) {
        TestIntegratedProbability test = new TestIntegratedProbability();
        NcmSimulator sim = new NcmSimulator();
        Network net = test.getScenario(4, 1);
        try {
            List<Network> nets = sim.createNetsWithNcmPrior(net, 100);
            List<Double> lengths = new ArrayList<>();
            for (Network n: nets) {
                lengths.addAll(sim.getBranchLengths(n));
            }
            Map<Double, Integer> histogram = sim.createBlHistogram(lengths, 50);
            System.out.println(histogram);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private List<Network> createNetsWithNcmPrior(Network topology, int num_settings) throws IOException {
        List<Network> nets = new ArrayList<>();
        String netFileString = "./tempNets.txt";

        Utils._START_NET = "[0.1]" + topology.toString();
        Utils._SAMPLE_FREQUENCY = 500;
        Utils._BURNIN_LEN = 1000 * Utils._SAMPLE_FREQUENCY;
        Utils._CHAIN_LEN = num_settings * Utils._SAMPLE_FREQUENCY + Utils._BURNIN_LEN;
        Utils._BL_EXP_PRIOR = true;
        Utils._TIMES_EXP_PRIOR = false;
        Utils.EXP_PARAM = 2 * _prior_bl_mean; // todo: figure out why this needs to be doubled?
        Utils.EXP_CUTOFF = _prior_bl_max;
        Utils._CONST_POP_SIZE = false;
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._ESTIMATE_POP_PARAM = false;
        Utils._NETWORK_SIZE_PRIOR = false;
        Utils._SEED = 1;
        Utils.DISABLE_TOPOLOGY_MOVES = true;

        MC3Core run = new MC3Core();

        // todo: use MC3Core._networkList, which is public (for now...?)
        PrintStream console = System.out;
        File netFile = new File(netFileString);
        FileOutputStream fos = new FileOutputStream(netFile);
        PrintStream ps = new PrintStream(fos);
        System.setOut(ps);
        // Run the MCMC network generation, writing stdout to treeFile
        run.run();
        // Return stdout to the console
        System.setOut(console);

        BufferedReader reader = new BufferedReader(new FileReader(netFile));
        String line = "";
        String netString;
        int lineNum = 0;
        int netNum;

        while (reader.ready() && ! line.contains("Summarization")) {
            line = reader.readLine();
            if (lineNum % 5 == 3) {
                netNum = (lineNum - 3) / 5;
            } else {netNum = -1;}

            if (netNum > Utils._BURNIN_LEN / Utils._SAMPLE_FREQUENCY) {
                netString = line.substring(line.indexOf("("));
                Network networkWithParams = Networks.readNetwork(netString);
                nets.add(networkWithParams);
//                System.out.println("net: " + netString);

//                System.out.println("produced trees: " + simulatedTrees.toString());
            }

            lineNum += 1;
        }
        return nets;
    }

    /**
     * Use MCMC to simulate networks with a given topology and distribution for branch length/inheritance probability
     * Then use those networks to simulate individual gene trees, matching the No Common Mechanism model
     * @param topology
     * @param num_settings counts the intermediate simulated networks (with same topology)
     * @return
     * @throws IOException
     */
    public List<Tree> simulateTrees(Network topology, int num_settings, int trees_per_setting) throws IOException {
        List<Tree> trees = new ArrayList<>();
        SimGTInNetworkByMS simulator = new SimGTInNetworkByMS();
        List<Network> nets = createNetsWithNcmPrior(topology, num_settings);
        for (Network networkWithParams: nets) {
            List<Tree> simulatedTrees = simulator.generateGTs(networkWithParams, null, trees_per_setting, "/Users/hunter/rice_grad/ms_folder/msdir/ms");
            trees.addAll(simulatedTrees);
        }
//        System.out.println("average bl was: " + average_bl);
//        _prior_bl_mean = average_bl;


        return trees;
    }

    private Map<Double, Integer> createBlHistogram(List<Double> bls, int numBins) {
        Map<Double, Integer> histogram = new HashMap<>();
        bls.sort(Comparator.<Double>naturalOrder());
        Double[] range = {bls.get(0), bls.get(bls.size() - 1)};
//        System.out.println(range[1]);
//        System.out.println(bls);
        Double[] bins = new Double[numBins];
        for (int i = 0; i < numBins; i++) {
            bins[i] = range[0] + i * (range[1] - range[0]) / numBins;
        }

        for (Double bl: bls) {
            int binId = 0;
//            System.out.printf("bin: %f, bl: %f", bins[binId], bl);
            while (binId < bins.length && bins[binId] <= bl) binId += 1;
            binId -= 1;
//            histogram.put(bins[binId], 1);
            histogram.put(bins[binId], 1 + histogram.getOrDefault(bins[binId], 0));
        }
        return histogram;
    }

    private List<Double> getBranchLengths(Network<Double> net) {
        List<Double> branch_lengths = new ArrayList<>();
        for (NetNode<Double> node: Networks.postTraversal(net)) {
//            if (node.isNetworkNode()) {
            node.getParents().forEach((NetNode parent) -> {

                branch_lengths.add(node.getParentDistance(parent));
            });
//                total += node.getParentDistance(node.getParents().)
//            }
        }
        return branch_lengths;
    }

    public Map<Tree, Double> getTopologyFrequencies(List<Tree> trees) {
        Map<Tree, Double> topologyCounts = new HashMap<>();
        for (Tree tree: trees) {
            boolean found = false;
            for (Tree match: topologyCounts.keySet()) {
                if (Trees.haveSameRootedTopology(tree, match)) {
                    topologyCounts.put(match, topologyCounts.get(match) + 1.0 / trees.size());
                    found = true;
                    break;
                }
            }
            if (! found) {
                // strip the tree's branch lengths
                for (TNode node: tree.getNodes()) {
                    if (node != tree.getRoot()) {
                        node.setParentDistance(TNode.NO_DISTANCE);
                    }
                }
                topologyCounts.put(tree, 1.0 / trees.size());
            }
        }
        return topologyCounts;
    }

    private double get_avg_branch_length(Network<Double> net) {
        double total = 0.0;
        for (Double bl: getBranchLengths(net)) {
            total += bl;
        }
        return total / net.getEdgeCount();
    }
}
