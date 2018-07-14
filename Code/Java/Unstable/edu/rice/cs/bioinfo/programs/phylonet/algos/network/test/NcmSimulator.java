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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by hunter on 7/12/18.
 */
public class NcmSimulator {
    public double _prior_bl_mean = 1.0;
    public double _prior_bl_max = Double.POSITIVE_INFINITY;
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
        String treeFileString = "./tempTrees.txt";

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

        // todo: is there a way to get results from MC3Core.run without this hack?
        // yes: use MC3Core._networkList, which is public (for now...?)
        PrintStream console = System.out;
        File treeFile = new File(treeFileString);
        FileOutputStream fos = new FileOutputStream(treeFile);
        PrintStream ps = new PrintStream(fos);
        System.setOut(ps);
        // Run the MCMC network generation, writing stdout to treeFile
        run.run();
        // Return stdout to the console
        System.setOut(console);

        BufferedReader reader = new BufferedReader(new FileReader(treeFile));
        String line = "";
        String netString;
        int lineNum = 0;
        int netNum;

        double average_bl = 0.0;

        SimGTInNetworkByMS simulator = new SimGTInNetworkByMS();

        while (reader.ready() && ! line.contains("Summarization")) {
            line = reader.readLine();
            if (lineNum % 5 == 3) {
                netNum = (lineNum - 3) / 5;
            } else {netNum = -1;}

            if (netNum > Utils._BURNIN_LEN / Utils._SAMPLE_FREQUENCY) {
                netString = line.substring(line.indexOf("("));
                Network networkWithParams = Networks.readNetwork(netString);

                double bl_contribution = get_avg_branch_length(networkWithParams);
                average_bl += bl_contribution;

//                System.out.println("net: " + netString + ", bl: " + bl_contribution);
                List<Tree> simulatedTrees = simulator.generateGTs(networkWithParams, null, trees_per_setting, "/Users/hunter/rice_grad/ms_folder/msdir/ms");
//                trees.add(simulatedTrees.get(0));
                trees.addAll(simulatedTrees);
//                System.out.println("produced trees: " + simulatedTrees.toString());
            }

            lineNum += 1;
        }

        average_bl /= num_settings;
//        System.out.println("average bl was: " + average_bl);
//        _prior_bl_mean = average_bl;


        return trees;
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
        List<Double> branch_lengths = new ArrayList<>();
        final double[] total = {0.0};
        for (NetNode<Double> node: Networks.postTraversal(net)) {
//            if (node.isNetworkNode()) {
            node.getParents().forEach((NetNode parent) -> {

                total[0] += node.getParentDistance(parent);
                branch_lengths.add(node.getParentDistance(parent));
            });
//                total += node.getParentDistance(node.getParents().)
//            }
        }
        return total[0] / net.getEdgeCount();
    }
}
