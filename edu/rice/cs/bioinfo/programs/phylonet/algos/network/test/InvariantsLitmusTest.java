package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityIntegrated;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import sun.nio.ch.Net;

import java.util.*;

/**
 * Created by hunter on 9/5/18.
 */
public class InvariantsLitmusTest {

    private static Comparator<List<Integer>> sortListByFirstElement = new Comparator<List<Integer>>() {
        public int compare(List<Integer> l1, List<Integer> l2) {
            return l1.get(0) - l2.get(0);
        }
    };

    public static void main(String[] args) {

        Network<Double> treeForAbba = Networks.readNetwork("(O, (G, (H, C)));");
        Network treeForAbbaWithoutOG = Networks.readNetwork("(C, (B, A));");
//        Network netForAbba = Networks.readNetwork();
        //abba baba tree case
        runTest(treeForAbba, "O");

        List<Tree> testTrees = getAllTreesMinusOutgroup(treeForAbba, "O");
        testTrees = getAllTreesWithOutgroup(treeForAbba, "O");

        //abba baba network case <- the values do not look right, needs to be looked into why
        //Network<Double> netForAbba = Networks.readNetwork("(O:1.0, ((H:1.0,(C:1.0)#H1:1::0.3):1.0,G:1.0,#H1:1::0.7));");
        double inheritanceProbability = 0.1;
        Network<Double> netForAbba = Networks.readNetwork("(O, ((G, #H1:0::" + inheritanceProbability + "), (H, (C)#H1:0::" + (1 - inheritanceProbability) + ")));");
        runTest(netForAbba,"O");

        //dfoil fig tree case
        Network<Double> treeForDFOIL = Networks.readNetwork("(((P1, P2), (P3, P4)), O);");
        //dfoil fig tree case
        runTest(treeForDFOIL, "O");


        //dfoil fig network case
        //add reticultion between P1 and P3 or P1 and P4 (can use whatever way hunter has been using on his own)
        inheritanceProbability = 0.1;
        Network<Double> netForDFOIL = Networks.readNetwork("(((P1, (P2)#H1:0::" + inheritanceProbability + "), ((P3, #H1:0::" + (1 - inheritanceProbability) + "), P4)), O);");
        //dfoil fig net case
        runTest(netForDFOIL, "O");



        //leo note - not really going past here for now

//        Network<Double> netForFullTest = Networks.readNetwork("(A, (B, (C, (D, E))));");
//        for (Tree t: getAllTreesMinusOutgroup(netForFullTest,"A")) {
//            Network netForTest = Networks.readNetwork(t.toNewick());
//            System.out.println();
//            runTest(netForTest,"A");
//        }
//        Network<Double> balancedNet = Networks.readNetwork("((A,C)hnode2,(B,D)hnode1);");
//        System.out.println(getAllTrees(balancedNet));
//        runTest(balancedNet);
    }

    /**
     *
     * @param topology is the tree or network topology (with inheritance probability), which contains outgroup
     * @param outgroupName
     */
    public static void runTest(Network topology, String outgroupName) {

        List<Tree> trees = getAllTreesWithOutgroup(topology, outgroupName);

        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
//        GeneTreeProbability gtp = new GeneTreeProbability();

        GeneTreeProbabilityIntegrated gtpi = new GeneTreeProbabilityIntegrated();
        double blParam = 1.0;
        double blMax = Double.POSITIVE_INFINITY;
        gtpi.setBranchLengthExponentialPrior(blParam, blMax);

        List<Network> parameterSettings = getAllParametrizedNetworks(topology, 1); // 100 random settings
        boolean firstSetting = true;
        // Invariant groups which are consistent across parameter settings
        List<List<Integer>> consistentInvariantGroups = new ArrayList<>();

        List<Double> probabilities = new ArrayList<>();
        for (Network netWithParameters: parameterSettings) {
            // use this line when gtp is a GeneTreeProbability
//            probabilities = gtp.calculateGTDistribution(netWithParameters, trees, null, false);
            // use these [6] lines when gtp is a GeneTreeProbabilityYF
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

        // Strip internal node names for convenience
        for (Object onode: topology.getTreeNodes()) {
            NetNode node = (NetNode)onode;
            if (!node.isLeaf()) node.setName("");
        }
        for (Object onode: topology.getNetworkNodes()) {
            NetNode node = (NetNode)onode;
            if (!node.isLeaf()) node.setName("");
        }

        String resultsString = "\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";
        resultsString += "\ntesting: " + topology.toString();
        resultsString += "\ntrees:";
        int i = 0;
        for (Tree t: trees) {
            resultsString += " " + i++ + ": " + t;
        }
        resultsString += "\ninvariants from random sampling and standard likelihood:";
        consistentInvariantGroups.sort(sortListByFirstElement);
        resultsString += "\n" + consistentInvariantGroups;
        resultsString += "\ninvariants from integrated likelihood:";
        resultsString += "\n" + integratedInvariantGroups;

        resultsString += "\n\nexample parameter setting:\n";
        resultsString += parameterSettings.get(parameterSettings.size() - 1);

        resultsString += "\nexample classic probs:\n";
        i = 0;
        for (Double d: probabilities) {
            resultsString += i++ + ": " + d + ", ";
        }

        resultsString += "\nintegrated probs:\n";
        i = 0;
        for (Double d: integratedProbs) {
            resultsString += i++ + ": " + d + ", ";
        }

        System.out.println(resultsString);
    }

    /**
     *
     * @param topo Topology which must include outgroupName as a leaf
     * @param outgroupName
     * @return A list of trees, which do not include the outgroup
     */
    private static List<Tree> getAllTreesMinusOutgroup(Network<Double> topo, String outgroupName) {
        // If this method throws an index out of bounds exception, ensure that topo contains outgroupName

        String[] theTaxa = new String[topo.getLeafCount()-1]; //{"A", "B", "C", "D"};
        int i = 0;
        for (NetNode leaf: topo.getLeaves()) {
            if(!leaf.getName().equals(outgroupName)){
                theTaxa[i++] = leaf.getName();
            }
        }
        Arrays.sort(theTaxa);

        List<Tree> oldTrees = Trees.generateAllBinaryTrees(theTaxa);
        List<Tree> trees = new ArrayList<>();

        for (Tree t: oldTrees) {
            // Recreate all trees so that node IDs will start at 0
            // Otherwise something goes wrong in GeneTreeProbability
            trees.add(Trees.readTree(t.toNewick()));
        }

        return trees;
    }

    /**
     *
     * @param topo Topology which must include outgroupName as a leaf
     * @param outgroupName
     * @return A list of trees, which include the outgroup
     */
    private static List<Tree> getAllTreesWithOutgroup(Network<Double> topo, String outgroupName) {
        List<Tree> withoutOutgroup = getAllTreesMinusOutgroup(topo, outgroupName);
        List<Tree> withOutgroup = new ArrayList<>();
        for (Tree t: withoutOutgroup) {
            withOutgroup.add(addOutgroup(t, outgroupName));
        }
        return withOutgroup;
    }

    // Check for fuzzy equality of (double) floating-point values
    private static boolean closeEnough(double x, double y) {
        double epsilon = 0.0001;
        double differenceMagnitude = Math.abs(x - y);
        double averageMagnitude = (Math.abs(x) + Math.abs(y)) / 2.0;
        return differenceMagnitude < epsilon * averageMagnitude;
    }

    // Find indices of (trees with) equal probabilities
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

    private static List<Network> getAllParametrizedNetworks(Network topology, int numRandomSettings) {
        List<Network> parametrizedNetworks = new ArrayList<>();

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
            for (NetNode parent: node.getParents()) {
                double newLength;
                newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//                newLength = 4.0;
//                newLength = Math.random() * 4; //10;
                node.setParentDistance(parent, newLength);

            }
        }
        for (NetNode<Double> hybridNode: newNet.getNetworkNodes()) {
            for (NetNode parent: hybridNode.getParents()) {
                double newLength;
                newLength = 0.0;
//                newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//                newLength = NetNode.NO_DISTANCE; // THIS IS NEGATIVE INFINITY, DO NOT USE
                hybridNode.setParentDistance(parent, newLength);
            }
        }
        return newNet;
    }

    private static Tree addOutgroup(Tree tree, String outgroupName) {
        TMutableNode oldRoot = (TMutableNode) tree.getRoot();
        STITree newTree = (STITree) Trees.readTree(";");
        STITree outgroup = (STITree) Trees.readTree(outgroupName);

        newTree.getRoot().adoptChild(oldRoot);
        newTree.getRoot().adoptChild(outgroup.getRoot());

        //TODO adopt child not working correctly in updating tree properties like set of nodes etc, adding this inefficient line to deal with it for now
        newTree = (STITree) Trees.readTree(newTree.toNewickWD());
        return newTree;
    }



    private static void FindGeneTreeInvariantsBracket(Tree treeHypothesis, Network networkHypothesis, String outgroupName){

    }

    private static void FindGeneTreeInvariantsIntegrated(Tree treeHypothesis, Network networkHypothesis, String outgroupName){

    }

    private static void TestDStatisticCase(){

    }

}