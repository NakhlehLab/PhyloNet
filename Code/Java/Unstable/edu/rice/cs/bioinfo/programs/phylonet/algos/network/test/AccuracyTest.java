package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferNetworkML;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferNetworkMLFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferNetworkNCM;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by hunter on 8/7/18.
 */
public class AccuracyTest {
    public static void main(String[] args) {
        TestIntegratedProbability test = new TestIntegratedProbability();
        AccuracyTest accTest = new AccuracyTest();
        List<Network> networksForTest = new ArrayList<>();
        for (int s = 1; s <= 4; s++) {
            for (double bl = 1; bl < 2.5; bl++) {
                networksForTest.add(test.getScenario(s, bl));
                networksForTest.add(test.getScenario(s, bl));
                networksForTest.add(test.getRandomNetwork());
            }
        }
        int[] dataSizes = {10, 100};
        int numTrials = 1;
        accTest.runInferences(networksForTest, dataSizes, numTrials);
    }

    public void runInferences(List<Network> networks, int[] dataSizes, int numTrials) {
        String msPath = "/Users/hunter/rice_grad/ms_folder/msdir/ms";
        SimGTInNetworkByMS simulator = new SimGTInNetworkByMS();
        InferNetworkNCM inferenceNcm = new InferNetworkNCM();
        InferNetworkML inferenceMl = new InferNetworkMLFromGTT_SingleTreePerLocus();
        for (int numGeneTrees: dataSizes) {
            System.out.printf("\n\ntrees: %d", numGeneTrees);

            int numSol = 10;
            double blParam = 1.0;
            int[] topN = new int[numSol];

            for (Network n: networks) {
                System.out.printf("\nnetwork: %s", n.toString());

                for (int trial = 0; trial < numTrials; trial++) {
                    int maxReticulations = ((List) n.getNetworkNodes()).size();
                    List<Tree> trees = simulator.generateGTs(n, null, numGeneTrees, msPath);
                    List<MutableTuple<Tree, Double>> dataNcm = getDataForNcmFromGts(trees);

                    List<Tuple<Network, Double>> results
                            = inferenceNcm.inferNetwork(dataNcm, null, maxReticulations, numSol, blParam);
                    //                System.out.printf("results: %s", results.toString());

                    for (int i = 0; i < numSol; i++) {
                        if (Networks.hasTheSameTopology(results.get(i).Item1, n)) {
                            topN[i] += 1;
                            System.out.printf(" placed: %d", i + 1);
                            break;
                        }
                    }
                }
            }

            System.out.printf("\nout of: %d", networks.size());
            int wins = 0;
            for (int i = 0; i < numSol; i++) {
                wins += topN[i];
                System.out.printf("\ntop %d: %d", i + 1, topN[i]);
            }
        }
    }

    private List<MutableTuple<Tree, Double>> getDataForNcmFromGts(List<Tree> trees) {
        List<MutableTuple<Tree, Double>> data = new ArrayList<>();
        for (Tree tree: trees) {
            boolean found = false;
            for (MutableTuple<Tree, Double> tuple: data) {
                Tree match = tuple.Item1;
                if (Trees.haveSameRootedTopology(tree, match)) {
                    tuple.Item2 += 1.0;
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
                data.add(new MutableTuple<Tree, Double>(tree, 1.0));
            }
        }
        return data;
    }
}
