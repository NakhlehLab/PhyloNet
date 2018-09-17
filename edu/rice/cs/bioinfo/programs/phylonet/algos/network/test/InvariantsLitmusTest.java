package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityIntegrated;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import sun.nio.ch.Net;

import java.util.*;

/**
 * Created by hunter on 9/5/18.
 */
public class InvariantsLitmusTest {

    public static void main(String[] args) {
        Network<Double> netWith5 = Networks.readNetwork("(A, (B, (C, (D, E))));");
        Network netWith6 = Networks.readNetwork("(A, (B, (C, (D, (E, F)))));");
//        runTest(netWith6);
        Network treeForAbba = Networks.readNetwork("(D, (C, (B, A)));");
//        Network netForAbba = Networks.readNetwork();
        runTest(treeForAbba);

        Network<Double> netForTest2 = Networks.readNetwork("(A, (B, (D, C)));");
//        System.out.println(getAllTrees(netForTest2));
        for (Tree t: getAllTrees(netWith5 )) {
            Network netForTest = Networks.readNetwork(t.toNewick());
//            runTest(netForTest);
        }
        Network<Double> balancedNet = Networks.readNetwork("((A,C)hnode2,(B,D)hnode1);");
//        System.out.println(getAllTrees(balancedNet));
//        runTest(balancedNet);



//        InvariantsLitmusTest.runTest(netForTest2);
//        InvariantsLitmusTest.test();
    }

    public static void test() {
        TestIntegratedProbability tip = new TestIntegratedProbability();

        GeneTreeProbability gtp = new GeneTreeProbability();
        GeneTreeProbabilityYF gtpyf = new GeneTreeProbabilityYF();
        GeneTreeProbabilityIntegrated gtpi = new GeneTreeProbabilityIntegrated();
        double blParam = 1.0;
        double blMax = Double.POSITIVE_INFINITY;
        gtpi.setBranchLengthExponentialPrior(blParam, blMax);

        Network s1bl1 = tip.getScenario(1, 1);
        Network net = randomizeBranchLengths(s1bl1);
        System.out.println("net: " + net);
        List<Tree> trees = getAllTrees(net);
        System.out.println("trees: " + trees);
        List<Double> gtpDist = gtp.calculateGTDistribution(net, trees, null, false);
        System.out.println("gtp: \t" + gtpDist);
        double[] gtpyfDist = new double[trees.size()];
        gtpyf.calculateGTDistribution(net, trees, null, gtpyfDist);
        System.out.println("gtpyf: \t" + Arrays.toString(gtpyfDist));
        List<Double> gtpiDist = gtpi.calculateGTDistribution(net, trees, null, false);
        System.out.println("gtpi: \t" + gtpiDist);

        List<Double> gtpyfDistList = new ArrayList<>();
        for (double d: gtpyfDist) {
            gtpyfDistList.add(d);
        }
    }


    public static void runTest(Network topology) {
        List<Tree> trees = getAllTrees(topology);

        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
//        GeneTreeProbability gtp = new GeneTreeProbability();

        GeneTreeProbabilityIntegrated gtpi = new GeneTreeProbabilityIntegrated();
        double blParam = 1.0;
        double blMax = Double.POSITIVE_INFINITY;
        gtpi.setBranchLengthExponentialPrior(blParam, blMax);

        List<Network> parameterSettings = getAllParametrizedNetworks(topology);
        boolean firstSetting = true;
//        List<Tuple<Integer, Integer>> consistentPairInvariants = new ArrayList<>();
        List<List<Integer>> consistentInvariantGroups = new ArrayList<>();

        List<Double> probabilities = new ArrayList<>();
        for (Network netWithParameters: parameterSettings) {
//            probabilities = gtp.calculateGTDistribution(netWithParameters, trees, null, false);
            double[] probsArray = new double[trees.size()];
            gtp.calculateGTDistribution(netWithParameters, trees, null, probsArray);
            probabilities = new ArrayList<>();
            for (double d: probsArray) {
                probabilities.add(d);
            }

            if (firstSetting) {
                firstSetting = false;
                consistentInvariantGroups = findInvariantGroups(probabilities);
            } else {
                consistentInvariantGroups = updateInvariantGroups(consistentInvariantGroups, probabilities);
            }
//            System.out.println(probabilities);
        }

        List<Double> integratedProbs = gtpi.calculateGTDistribution(topology, trees, null, false);

        List<List<Integer>> integratedInvariantGroups = findInvariantGroups(integratedProbs);

//        Networks.removeAllParameters(topology);
        String netString = topology.toString();
        String resultsString = "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";
        resultsString += "\ntesting: " + netString;
        resultsString += "\ntrees:";
        int i = 0;
        for (Tree t: trees) {
            resultsString += " " + i++ + ": " + t;
        }
        resultsString += "\ninvariants from random sampling and standard likelihood:";
        resultsString += "\n" + consistentInvariantGroups;
        resultsString += "\ninvariants from integrated likelihood:";
        resultsString += "\n" + integratedInvariantGroups;
//        resultsString += "\n\nclassic probs:\n";
//        int i = 0;
//        for (Double d: probabilities) {
//            resultsString += i++ + ": " + d + ", ";
//        }
//        resultsString += "\nintegrated probs:\n";
//        i = 0;
//        for (Double d: integratedProbs) {
//            resultsString += i++ + ": " + d + ", ";
//        }
        System.out.println(resultsString);
    }

    private static List<Tree> getAllTrees(Network<Double> topo) {

        String[] theTaxa = new String[topo.getLeafCount()]; //{"A", "B", "C", "D"};
        int i = 0;
        for (NetNode leaf: topo.getLeaves()) {
            theTaxa[i++] = leaf.getName();
        }

        List<Tree> oldTrees = Trees.generateAllBinaryTrees(theTaxa);
        List<Tree> trees = new ArrayList<>();

        for (Tree t: oldTrees) {
            // Recreate all trees so that node IDs will start at 0
            // Otherwise something goes wrong in the probability calculation
            trees.add(Trees.readTree(t.toNewick()));
        }

        return trees;
    }

    private static boolean closeEnough(double x, double y) {
        double epsilon = 0.0001;
        double differenceMagnitude = Math.abs(x - y);
        double averageMagnitude = (Math.abs(x) + Math.abs(y)) / 2.0;
        return differenceMagnitude < epsilon * averageMagnitude;
    }

    private static List<List<Integer>> findInvariantGroups(List<Double> probabilities) {
        List<List<Integer>> groups = new ArrayList<>();

        for (int i = 0; i < probabilities.size(); i++) {
            boolean foundAGroup = false;
            for (List<Integer> group: groups) {
                if (closeEnough(probabilities.get(i), probabilities.get(group.get(0)))) {
                    group.add(i);
                    foundAGroup = true;
                    break;
                }
            }
            if (! foundAGroup) {
                List<Integer> newGroup = new ArrayList<>();
                newGroup.add(i);
                groups.add(newGroup);
            }
        }

        return groups;
    }

    private static List<List<Integer>> updateInvariantGroups(List<List<Integer>> oldGroups, List<Double> probabilities) {
        List<List<Integer>> newGroups = new ArrayList<>();
        for (List<Integer> group: oldGroups) {
            while (group.size() > 0) {
                List<Integer> newGroup = new ArrayList<>();
                newGroup.add(group.get(0));
                group.remove(0);
                for (Integer i: group) {
                    if (closeEnough(probabilities.get(i), probabilities.get(newGroup.get(0)))) {
                        newGroup.add(i);
                    }
                }
                for (Integer i: newGroup) {
                    if (group.contains(i)) {
                        group.remove(group.indexOf(i));
                    }
                }
                newGroups.add(newGroup);
            }
        }
        return newGroups;
    }

    private static List<Network> getAllParametrizedNetworks(Network topology) {
        List<Network> parametrizedNetworks = new ArrayList<>();

        int numRandomSettings = 100;
        for (int i = 0; i < numRandomSettings; i++) {
            Network<Double> newNet = randomizeBranchLengths(topology);
            parametrizedNetworks.add(newNet);
        }
        return parametrizedNetworks;
    }

    private static Network randomizeBranchLengths(Network topology) {
        double[] branchLengths = {0.5, 1, 2, 4};
        Network<Double> newNet = topology.clone();

        for (NetNode<Double> node: newNet.getTreeNodes()) {
            double newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//            newLength = 4.0;
//            double newLength = Math.random() * 4; //10;
            for (NetNode parent: node.getParents()) {
                node.setParentDistance(parent, newLength);
            }
        }
        for (NetNode<Double> hybridNode: newNet.getNetworkNodes()) {

            for (NetNode parent: hybridNode.getParents()) {
                double newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//                newLength = 0.0;
//                newLength = NetNode.NO_DISTANCE; // THIS IS NEGATIVE INFINITY, DO NOT USE

                hybridNode.setParentDistance(parent, newLength);
            }
        }
        return newNet;
    }

    // Refactor these two for more DRY
    private static Network addOutgroup(Network net) {
        String outgroupName = "OG";
        // todo add outgroup
        return net;
    }

    // Refactor these two for more DRY
    private static Tree addOutgroup(Tree tree) {
        String outgroupName = "OG";
        // todo add outgroup
        return tree;
    }

}
