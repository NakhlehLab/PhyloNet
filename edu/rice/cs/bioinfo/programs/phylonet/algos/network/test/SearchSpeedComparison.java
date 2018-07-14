package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.*;

/**
 * Created by hunter on 7/5/18.
 *
 * run the same search using ML and ML-NCM hill climbing; compare how many steps were taken
 *
 */
public class SearchSpeedComparison {
    TestIntegratedProbability test = new TestIntegratedProbability();
    NcmSimulator ncmSim = new NcmSimulator();

    List<MutableTuple<Tree,Double>> inputData;
    Network trueNetwork;
    List<Tuple<Network,Double>> inferenceResultsNCM;
    LinkedList<Tuple<Network,Double>> inferenceResultsML;

    InferNetworkNCM inferNCM;
    InferNetworkML inferML;
    File NCMlogFile = new File("ncmlog.txt");
    File MLlogFile = new File("mllog.txt");

    public static void main(String[] args) throws IOException {
        SearchSpeedComparison search = new SearchSpeedComparison();
//        search.runSearches();
        for (int s = 1; s <= 4; s++) {
            for (int b = 1; b <= 2; b++) {
                int hybrids = 1;
                if (s > 1) hybrids += 1;
                System.out.println("\n\n\n////////////////////////////////////\ns" + s + "bl" + b + "(" + hybrids + "):\n////////////\n");
                search.runSearches(s, b, hybrids);
            }
        }
    }

    public void setData(List<MutableTuple<Tree,Double>> data) {
        inputData = data;
    }

    public void runSearches(int scenN, double bl, int hybrids) throws IOException {
        double blMean = bl;
        Network<Double> trueNetwork = test.getScenario(scenN, blMean);
        int numReticulations = hybrids;
        int numResults = 2;
        ncmSim._prior_bl_max = Double.POSITIVE_INFINITY;
        ncmSim._prior_bl_mean = blMean;

        List trees = ncmSim.simulateTrees(trueNetwork, 100, 1);
        Map<Tree, Double> treeFrequencies = ncmSim.getTopologyFrequencies(trees);
        List<MutableTuple<Tree, Double>> data = new ArrayList<>();
        for (Tree topology: treeFrequencies.keySet()) {
            data.add(new MutableTuple<>(topology, treeFrequencies.get(topology)));
        }

        setData(data);

        // search parameters

        int maxRounds = 100;
        int maxTryPerBranch = 100;
        double improvementThreshold = 0.001;
        double maxBranchLength = 6;
        double Brent1 = 0.01;
        double Brent2 = 0.005;
        long maxExaminations = 1000000;
        int maxFailure = 100;
        int moveDiameter = -1;
        int reticulationDiameter = -1;
        int numProcessors = 1;
        Network startNetwork = null;
        Set<String> fixedHybrid = new HashSet<String>();;
        double[] topologyOperationWeights = {0.1,0.1,0.15,0.55,0.15,0.15};
        double[] topologyVsParameterOperation = {0.3,0.7}; // modify operationWeights if you want pure parameter operation
        double[] operationWeights = {0.1,0.1,0.15,0.55,0.15,0.15, topologyVsParameterOperation[1] / topologyVsParameterOperation[0]}; // don't devide by 0
        int numRuns = 2;
        boolean optimizeBL = false;
        Long seed = Long.valueOf(100);



        System.out.println("\n\nbeginning ncm inference");

        inferNCM = new InferNetworkNCM();
//        inferNCM.setLogFile(NCMlogFile);
        inferNCM.setSearchParameter(maxExaminations, maxFailure, moveDiameter, reticulationDiameter, startNetwork, fixedHybrid, numProcessors, topologyOperationWeights, numRuns, seed);
        inferenceResultsNCM = inferNCM.inferNetwork(inputData, null, numReticulations, numResults, blMean);

        List<List<MutableTuple<Tree, Double>>> dataForML = new ArrayList<>();
        dataForML.add(inputData);
        inferenceResultsML = new LinkedList<>();

        System.out.println("\n\nbeginning ml inference");

        inferML = new InferNetworkMLFromGTT_SingleTreePerLocus();
//        inferML.setLogFile(MLlogFile);
        inferML.setSearchParameter(maxRounds, maxTryPerBranch, improvementThreshold, maxBranchLength, Brent1, Brent2, maxExaminations, maxFailure, moveDiameter, reticulationDiameter, numProcessors, startNetwork, fixedHybrid, operationWeights, numRuns, optimizeBL, seed);
        inferML.inferNetwork(dataForML, null, numReticulations, numResults, false, inferenceResultsML);

//        System.out.println("printing results");
//
//        printResults();
    }

    public void printResults() {
//        System.out.println("data: " + inputData);
        System.out.println("ncm results: " + inferenceResultsNCM);
        System.out.println("ml results: " + inferenceResultsML);
    }
}
