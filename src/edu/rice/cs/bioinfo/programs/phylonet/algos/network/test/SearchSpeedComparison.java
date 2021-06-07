package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.nio.Buffer;
import java.util.*;

/**
 * Created by hunter on 7/5/18.
 *
 * run the same search using ML and ML-NCM hill climbing; compare the search time
 *
 */
public class SearchSpeedComparison {
    static TestIntegratedProbability test = new TestIntegratedProbability();
    static NcmSimulator ncmSim = new NcmSimulator();

    List<MutableTuple<Tree,Double>> inputData;
    Network trueNetwork;
    List<Tuple<Network,Double>> inferenceResultsNCM;
    LinkedList<Tuple<Network,Double>> inferenceResultsML;

    InferNetworkNCM inferNCM;
    InferNetworkML inferML;
    File NCMlogFile = new File("ncmlog.txt");
    File MLlogFile = new File("mllog.txt");

//    static String msPath = "edu/rice/cs/bioinfo/programs/phylonet/algos/network/test/ms";
    static String msPath = "./ms";

    public static void main(String[] args) throws IOException {
        int gtCount = 1000;
        SearchSpeedComparison search = new SearchSpeedComparison();
//        search.runSearches();
//        for (int scen = 1; scen <= 4; scen++) {
//            for (int bl = 1; bl <= 2; bl++) {
        for (int scen = 1; scen <= 4; scen++) {
            for (int bl = 1; bl <= 2; bl++) {
                Network trueNetwork = test.getScenario(scen, bl);
                int hybrids = trueNetwork.getReticulationCount();
                System.out.println("\n\n\n////////////////////////////////////\ns" + scen + "bl" + bl + "(h" + hybrids + "):\n////////////\n");

                ///////// Simulation
//                ncmSim._prior_bl_max = Double.POSITIVE_INFINITY;
//                ncmSim._prior_bl_mean = bl;
//                List trees = ncmSim.simulateTrees(trueNetwork, 100, 1);
                SimGTInNetworkByMS simulator = new SimGTInNetworkByMS();
                List trees = simulator.generateGTs(trueNetwork, null, gtCount, msPath);
                ///////// End Simulation

                Map<Tree, Double> treeFrequencies = ncmSim.getTopologyFrequencies(trees);
                search.runSearches(trueNetwork, treeFrequencies, hybrids, bl);
            }
        }
//        search.runSearches(1, 1.0, 1);
    }

    public void setData(List<MutableTuple<Tree,Double>> data) {
        inputData = data;
    }

//    public void runSearches(int scenN, double bl, int hybrids) throws IOException {
    public void runSearches(Network<Double> trueNetwork, Map<Tree, Double> treeFrequencies, int hybrids, double bl) throws IOException {
        int numResults = 6;

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
        int maxFailure = 1000; // 100 is default
        int moveDiameter = -1;
        int reticulationDiameter = -1;
        int numProcessors = 1;
        Network startNetwork = null;
        Set<String> fixedHybrid = new HashSet<String>();;
        double[] topologyOperationWeights = {0.1,0.1,0.15,0.55,0.15,0.15};
        double[] topologyVsParameterOperation = {0.3,0.7}; // modify operationWeights if you want pure parameter operation
        double[] operationWeights = {0.1,0.1,0.15,0.55,0.15,0.15, topologyVsParameterOperation[1] / topologyVsParameterOperation[0]}; // don't divide by 0
        int numRuns = 1;
        boolean optimizeBL = false; // this sets search to hillclimb in infernetML
        Long seed = 300L;



        System.out.println("\n\nbeginning ncm inference");

        inferNCM = new InferNetworkNCM();
        clearLogFile(NCMlogFile);
        inferNCM.setLogFile(NCMlogFile);
        inferNCM.setSearchParameter(maxExaminations, maxFailure, moveDiameter, reticulationDiameter, startNetwork, fixedHybrid, numProcessors, topologyOperationWeights, numRuns, seed);
        inferenceResultsNCM = inferNCM.inferNetwork(inputData, null, hybrids, numResults, bl);

        checkValidity(inferenceResultsNCM, trueNetwork);
        System.out.println(countLines(NCMlogFile) + " lines (networks searched)");



        System.out.println("\n\nbeginning ml inference");

        List<List<MutableTuple<Tree, Double>>> dataForML = new ArrayList<>();
        dataForML.add(inputData);
        inferenceResultsML = new LinkedList<>();

        inferML = new InferNetworkMLFromGTT_SingleTreePerLocus();
        clearLogFile(MLlogFile);
        inferML.setLogFile(MLlogFile);
        inferML.setSearchParameter(maxRounds, maxTryPerBranch, improvementThreshold, maxBranchLength, Brent1, Brent2, maxExaminations, maxFailure, moveDiameter, reticulationDiameter, numProcessors, startNetwork, fixedHybrid, operationWeights, numRuns, optimizeBL, seed);
        inferML.inferNetwork(dataForML, null, hybrids, numResults, false, inferenceResultsML);

        checkValidity(inferenceResultsML, trueNetwork);
        System.out.println(countLines(MLlogFile) + " lines (networks searched)");
//        System.out.println("printing results");
//
//        printResults();
    }

    private int countLines(File f) {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(f));
            int lines = 0;
            while (reader.readLine() != null) {
                lines += 1;
            }
            return lines;

        } catch(Exception e) {
            System.err.print(e.getStackTrace());
            return -1;
        }

    }

    private void clearLogFile(File f) {
        if (f != null) {
            f.delete();
        }
    }

    private void checkValidity(List<Tuple<Network, Double>> results, Network trueNetwork) {
//        boolean somewhereOutThere = false;
//        for (Tuple<Network, Double> t : results) {
//            if (Networks.hasTheSameTopology(trueNetwork, t.Item1)) {
//                somewhereOutThere = true;
//            }
//        }
//        boolean noPointsForSecondPlace = Networks.hasTheSameTopology(results.get(0).Item1, trueNetwork);
//        System.out.println("Correct network found: " + somewhereOutThere + "\nCorrect network first: " + noPointsForSecondPlace);
        String placement = ">= " + results.size();
        for (int i = 0; i < results.size(); i++) {
            Tuple<Network, Double> t = results.get(i);
            if (Networks.hasTheSameTopology(trueNetwork, t.Item1)) {
                placement = "" + i;
                break;
            }
        }
        System.out.println("Correct network found at position: " + placement);
    }

    public void printResults() {
//        System.out.println("data: " + inputData);
        System.out.println("ncm results: " + inferenceResultsNCM);
        System.out.println("ml results: " + inferenceResultsML);
    }
}
