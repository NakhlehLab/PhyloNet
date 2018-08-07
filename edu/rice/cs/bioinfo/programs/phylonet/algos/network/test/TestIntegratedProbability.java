package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.core.MC3Core;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTreeEnumerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.*;


/**
 * Created by hunter on 4/4/18.
 */
public class TestIntegratedProbability {

    private double _prior_bl_mean = 1;
    private double _prior_bl_max = Double.POSITIVE_INFINITY; // don't change this until you've updated GTPIntegrated

    private static GeneTreeProbabilityIntegrated gtpi = new GeneTreeProbabilityIntegrated();

    public static void main(String[] args) throws IOException {
        TestIntegratedProbability test = new TestIntegratedProbability();
        SearchSpeedComparison compare = new SearchSpeedComparison();

//        compare.runSearches();
//        test.troubleshoot();
        String[] stringsToClean = {};
        for (String toClean: stringsToClean) {
            Network n = Networks.readNetwork(toClean);
            Networks.removeAllParameters(n);
            System.out.println(n);
        }
//        Iterable<NetworkTree> enumerator = new NetworkTreeEnumerator(test.getScenario(4, 1));
        Iterable<NetworkTree> altIterable = Networks.getTrees(test.getScenario(4, 1));
//        Iterator<Tree> iter = enumerator.iterator();
        for (NetworkTree nt: altIterable) {
            System.out.println(nt.makeTree());
        }
    }

    private static void troubleshoot() {
        Map a2s = null;
        boolean printDetails = true;
        Network<Double> net = Networks.readNetwork("(A,((B)I3#H1:1.0::0.5,((I3#H1:::0.5,C)I1,D)I2)I4)I0;");
        List<Tree> gts = new ArrayList();
        gts.add(Trees.readTree("(A, (B, (C, D)));"));
        gts.add(Trees.readTree("(A, ((B, C), D));"));
        gtpi.calculateGTDistribution(net, gts, a2s, printDetails);
        net = Networks.readNetwork("(A,(((C,(D)I4#H1:1.0::0.5)I3,B)I2,I4#H1:::0.5)I1)I0;");
        gtpi.calculateGTDistribution(net, gts, a2s, printDetails);
        net = Networks.readNetwork("(A,(((C)I3#H1:::0.5,B)I2,(I3#H1:1.0::0.5,D)I4)I1)I0;");
        gtpi.calculateGTDistribution(net, gts, a2s, printDetails);
    }

    private static boolean isUltrametric(Network n) {
        double eps = 0.0001;
        Iterable<NetNode> leaves = n.getLeaves();
        Double netHeight = null;
        for (NetNode leaf: leaves) {
            if (netHeight == null) {
                netHeight = getHeight(n, leaf);
            } else if (Math.abs(netHeight - getHeight(n, leaf)) > eps) {
                return false;
            }
        }
        return true;
    }

    private static double getHeight(Network n, NetNode node) {
        if (node == n.getRoot()) {
            return 0;
        } else {
            Iterable<NetNode> parents = node.getParents();
            Double height = null;
            for (NetNode p: parents) {
                if (height == null) {
                    height = node.getParentDistance(p) + getHeight(n, p);
                } else if (height != node.getParentDistance(p) + getHeight(n, p)) {
                    return Double.NEGATIVE_INFINITY; // indicates a non-ultrametric network
                }
            }
            return height;
        }
    }

    /**
     * Networks for test as seen in section 3.4 of Yun Yu's thesis
     * @param i
     * @param bl
     * @return
     */
    public Network getScenario(int i, double bl) {
        String networkString = "";
        Network scenarioNetwork;
        double t1, t2, t3, t4, t5, t6;
        double ih1, ih2;
        t1 = bl;
        t2 = bl;
        t3 = bl;
        t4 = bl;
        t5 = bl;
        t6 = bl;
        ih1 = 0.5;
        ih2 = 0.5;
        switch (i) {
            case 1:
                networkString = "((A:" + (t1+t2+t3) + ", ((B:" + (t1) + ", C:" + (t1) + "):" + (t2) + ")#H1:" + (t3) + "::" + (ih1) + "):" + (t4) + ", (#H1:" + (t3) + "::" + (1 - ih1) + ", D:" + (t1+t2+t3) + "):" + (t4) + ");";
                break;
            case 2:
                String[] leaves = {"A", "B", "C", "D", "E", "F"};
                String[] branches = new String[2];
                double[] ihs = {ih1, ih2};
                for (int b = 0; b < 2; b++) {
                    branches[b] = "((" + leaves[0 + 3 * b] + ":" + (t1 + t2) + ", (" + leaves[1 + 3 * b] + ":" + t1 + ")#H" + (b + 1) + ":" + (t2) + "::" + (ihs[b]) + "):" + (t3) + ", (#H" + (b + 1) + ":" + (t2) + "::" + (1 - ihs[b]) + ", " + leaves[2 + 3 * b] + ":" + (t1 + t2) + "):" + t3 + "):" + t4;
                }
                networkString = "(" + branches[0] + ", "+ branches[1] + ");";
                break;
            case 3:
                networkString = "((A:" + (t1 + t2 + t3 + t4 + t5) + ", (((B:" + (t1 + t2) + ", (C:" + (t1) + ")#H2:" + (t2) + "::" + (ih2) + "):" + (t3) + ", (#H2:" + (t2) + "::" + (1 - ih2) + ", D:" + (t1 + t2) + "):" + (t3) + "):" + (t4) + ")#H1:" + (t5) + "::" + (ih1) + "):" + (t6) + ", (#H1:" + (t5) + "::" + (1 - ih1) + ", E:" + (t1 + t2 + t3 + t4 + t5) + "):" + (t6) + ");";
                break;
            case 4:
                String h2string = "(D:" + (t1) + ")";
                String h1string = "((B:" + (t1 + t2 + t3) + ", (C:" + (t1 + t2) + ", " + h2string + "#H2:" + (t2) + "::" + (ih2) + "):" + (t3) + "):" + (t4) + ")";
                String nodeAH = "(A:" + (t1 + t2 + t3 + t4 + t5) + ", " + h1string + "#H1:" + (t5) + "::" + (ih1) + "):" + (t6);
                String nodeEF = "((#H2:" + (t2) + "::" + (1 - ih2) + ", E:" + (t1 + t2) + "):" + (t3) + ", F:" + (t1 + t2 + t3) + "):" + (t4 + t5);
                String nodeHEF = "(#H1:" + (t5) + "::" + (1 - ih1) + ", " + nodeEF + "):" + (t6);
                networkString = "(" + nodeAH + ", " + nodeHEF + ");";
                break;
            default:
                throw new IllegalArgumentException();
        }
        scenarioNetwork = Networks.readNetwork(networkString);
        return scenarioNetwork;
    }

    private List<Tree> extractTrees(BufferedReader br, int treeCount) throws IOException {
        List<Tree> treeList = new ArrayList<>();
        for (int i=0; i < 4; i++) {br.readLine();} // trees start on 5th line
        for (int i=0; i < treeCount; i++) {
            String treeString = br.readLine();
            treeString = treeString.substring(treeString.indexOf('('), treeString.length());
            Tree t = Trees.readTree(treeString);
            treeList.add(t);
        }
        return treeList;
    }

    /**
     * read result networks from the output of an inference
     * @param br
     * @param networkCount
     * @return
     * @throws IOException
     */
    private List<Tuple<Network, Double>> extractResults(BufferedReader br, int networkCount) throws IOException {
        List output = new ArrayList();
        // 0 skipped lines before networks
        for (int i = 0; i < networkCount; i++) {
            br.readLine();
            br.readLine();
            // 2 skipped lines between networks
            Network n = Networks.readNetwork(br.readLine());
            String probLine = br.readLine();
            double prob = Double.parseDouble(probLine.substring("Total log probability:".length(), probLine.length()));
            output.add(new Tuple(n, prob));
        }
        return output;
    }

    /**
     *
     * @throws IOException
     */
    private void validateExperiment() throws IOException {
        String date_of_test = "6_11";
        String resultsLocation = "/Users/hunter/rice_grad/results_from_nots_" + date_of_test + "/";
        String dataLocation = "/Users/hunter/rice_grad/for_nots_" + date_of_test + "/";
        String testOutLocation = "/Users/hunter/rice_grad/validation_results_" + date_of_test + ".txt";
        String plotDataLocation = "/Users/hunter/rice_grad/plot_data_" + date_of_test + ".txt";
        int[] scenarios = {1, 2, 3, 4};
        double[] blSettings = {1, 2};
        int[] treeCounts = {10, 30, 100, 300, 1000, 3000};
        int trials = 20;
        NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
//        GeneTreeProbabilityIntegrated gtpi = new GeneTreeProbabilityIntegrated();

        Network trueNetwork;
        Network n;
        String probLine;
        double logProb = 0;
        double correctLogProb = 0;
        double averageDistance;
        double networkDistance;
        int perfectCount;
        int blameTheSearch;
        int trialsWithDistance;
        Boolean trueNetworkInTopN = false;
        int topN;
        int resultsLength = 0;

        // todo: write to file in a cleaner way
        PrintStream console = System.out;
        File outputFile = new File(testOutLocation);
        FileOutputStream fos = new FileOutputStream(outputFile);
        PrintStream ps = new PrintStream(fos);
        System.setOut(ps);

        // perhaps like this?
        FileWriter plotWriter = new FileWriter(plotDataLocation);
        plotWriter.write("{\"title\": \"Data for accuracy plot from 6/11\", \"data\": [");

        for (int s: scenarios) {
            for (double bl: blSettings) {
                trueNetwork = getScenario(s, bl);

                plotWriter.write("{\"scenario\": " + s + ", \"branch_length\": " + bl + ", \"data_for_plot\": [");

                for (int tc: treeCounts) {
                    averageDistance = 0;
                    perfectCount = 0;
                    blameTheSearch = 0;
                    trialsWithDistance = 0;
                    topN = 0;
                    boolean someNetworkBeatCorrect = false;

                    plotWriter.write("{\"tree_count\": " + tc + ", \"distances\": [");
                    String distancesString = "";


                    for (int t = 0; t < trials; t++) {



                        String resultsPath = resultsLocation + "s" + s + "trees" + tc + "bl" + (int)bl + "_" + t + ".out";
                        String dataPath = dataLocation + "scen" + s + "trees" + tc + "bl" + (int)bl + "_" + t + ".nex";
                        // try-with-resources
                        try (FileReader fr = new FileReader(resultsPath);
                                BufferedReader resultsReader = new BufferedReader(fr);
                                FileReader dataFileReader = new FileReader(dataPath);
                                BufferedReader dataReader = new BufferedReader(dataFileReader);) {

                            List<Tuple<Network, Double>> results = extractResults(resultsReader, 5);
                            resultsLength = results.size();
                            n = results.get(0).Item1;
                            logProb = results.get(0).Item2;
                            if (Networks.hasTheSameTopology(n, trueNetwork)) {
                                perfectCount += 1;
                                blameTheSearch += 1;
                                topN += 1;
                                distancesString += "0, ";
                            } else {
                                networkDistance = metric.computeDistanceBetweenTwoNetworks(trueNetwork, n);
                                distancesString += networkDistance + ", ";

                                trialsWithDistance += 1;
                                averageDistance += networkDistance;

                                List<Tree> treesForInference = extractTrees(dataReader, tc);
                                gtpi.setBranchLengthExponentialPrior(bl, Double.POSITIVE_INFINITY);
//                                List<Double> treeProbabilities = gtpi.calculateGTDistribution(trueNetwork, treesForInference, null, false);
//                                correctLogProb = 0.;
//                                for (Double prob : treeProbabilities) {
//                                    correctLogProb += Math.log(prob);
//                                }
                                correctLogProb = gtpi.getTotalLogLikelihood(trueNetwork, treesForInference);
                                if (correctLogProb < logProb) {
                                    //System.out.println("trial " + t + " failed with " + correctLogProb + " < " + logProb);
                                    trueNetworkInTopN = false;
                                    for (Tuple<Network, Double> tup: results) {
                                        if (Networks.hasTheSameTopology(tup.Item1, trueNetwork)) {
                                            trueNetworkInTopN = true;
                                        }
                                    }
                                    if (trueNetworkInTopN) topN += 1;
                                } else {
                                    blameTheSearch += 1;
                                }
                            }
                        }

                    }
                    if (trialsWithDistance > 0) {
                        averageDistance = averageDistance / trialsWithDistance;
                    }
                    System.out.println("s" + s + " bl" + bl + " with " + tc + " trees has " + perfectCount + " out of " + trials + " trials perfect, and average error (of incorrect nets): " + averageDistance);
                    System.out.println("but " + blameTheSearch + " trials give higher probability to the true network");
                    System.out.println("and " + topN + "trials found true network in the top " + resultsLength);
                    System.out.println();

                    // deal with trailing comma
//                    if (distancesString.contains(",")) {
//                        distancesString = distancesString.substring(0, distancesString.lastIndexOf(","));
//                    }
                    plotWriter.write(distancesString + "]}, ");
                }
                String settingDelim = "\n";
                for (int i = 0; i < 79; i++) settingDelim += "/";
                settingDelim += "\n";
                System.out.println(settingDelim);

                plotWriter.write("]}, ");
            }
        }

        plotWriter.write("]}");

        System.setOut(console);
        plotWriter.close();



    }

    private void testFrequencyOfML() throws IOException {
        int num_settings = 1;
        int trees_per_setting = 100000;
        Network tree1 = Networks.readNetwork("((A:4,B:4):1,C:5);");
        Network net2 = Networks.readNetwork("((A:2,(B:1)#H1:1::0.3):1,(#H1:1::0.7,C:2):1);");
        Network s4bl1 = Networks.readNetwork("((A:3, ((B:2, (C:1, (D:1)#H2:0::0.5):1):1)#H1:0::0.5):1, (#H1:0::0.5, ((#H2:0::0.5, E:1):1, F:2):1):1);");
        Network bigTree = Networks.readNetwork("(A:4, (B:3, (C:2, (D:1, E:1):1):1):1);");
        Network jtree = Networks.readNetwork("(((a:1, b:1):1, c:2):1, d:3);");
        Network[] networksForTest = {s4bl1};

        SimGTInNetworkByMS simulator = new SimGTInNetworkByMS();
        NcmSimulator ncmSim = new NcmSimulator();

        for (Network net: networksForTest) {

            List<Tree> simulatedTrees = simulator.generateGTs(net, null, trees_per_setting, "/Users/hunter/rice_grad/ms_folder/msdir/ms");

            Map<Tree, Double> frequencies = ncmSim.getTopologyFrequencies(simulatedTrees);

            GeneTreeProbability gtp = new GeneTreeProbability();

            List<Tree> simulatedTopologies = new ArrayList<>();
            for (Tree t : frequencies.keySet()) {
                simulatedTopologies.add(t);
            }

            double mean_error = 0.0;
            double reference_error = 0.0;

            List<Double> probabilities = gtp.calculateGTDistribution(net, simulatedTopologies, null, false);
            List<Double> pdf = normalizeDistribution(probabilities);
            List<Double> cdf = pdfToCdf(pdf);

            Map<Tree, Double> reference = new HashMap<>();
            int index;
            for (Tree top : simulatedTopologies) {
                reference.put(top, 0.0);
            }
            for (Tree st : simulatedTrees) {
                index = randomChoice(cdf);
                reference.put(simulatedTopologies.get(index), 1.0 / simulatedTrees.size() + reference.get(simulatedTopologies.get(index)));
            }

//            System.out.println(probabilities);
            int i = 0;
            for (Tree tree : simulatedTopologies) {
                System.out.println(tree.toString() + ", freq: " + frequencies.get(tree) + ", prob: " + probabilities.get(i) + ", ref: " + reference.get(tree));
                double ratio = frequencies.get(tree) / probabilities.get(i);
                if (ratio < 1) {
                    ratio = 1.0 / ratio;
                }
                mean_error += ratio;
                System.out.print("freq error got " + ratio);

                ratio = reference.get(tree) / probabilities.get(i++);
                if (ratio < 1) {
                    ratio = 1.0 / ratio;
                }
                if (ratio < Double.POSITIVE_INFINITY) {
                    reference_error += ratio;
                    System.out.println(", reference error got " + ratio);
                } else {
                    System.out.println(", tree " + tree.toString() + " did not appear in reference sample");
                }
            }
            mean_error /= simulatedTopologies.size();
            reference_error /= simulatedTopologies.size();
            System.out.println("mean error: " + mean_error);
            System.out.println("reference error: " + reference_error);
        }
    }

    private void testFrequencyAgainstProbability() throws IOException {
        int num_settings = 1000;
        int trees_per_setting = 1;
        Network tree1 = Networks.readNetwork("((A:4,B:4):1,C:5);");
        Network net2 = Networks.readNetwork("((A:2,(B:1)#H1:1::0.3):1,(#H1:1::0.7,C:2):1);");
        Network s4bl1 = Networks.readNetwork("((A:3, ((B:2, (C:1, (D:1)#H2:0::0.5):1):1)#H1:0::0.5):1, (#H1:0::0.5, ((#H2:0::0.5, E:1):1, F:2):1):1);");
        Network bigTree = Networks.readNetwork("(A:4, (B:3, (C:2, (D:1, E:1):1):1):1);");
        Network jtree = Networks.readNetwork("(((a:1, b:1):1, c:2):1, d:3);");
        Network s2bl1 = Networks.readNetwork("((E:5.0,(((D:2.0,(C:1.0)#H2:1.0::0.5):1.0,(#H2:1.0::0.5,B:2.0):1.0):1.0)#H1:1.0::0.5):1.0,(#H1:1.0::0.5,A:5.0):1.0);");
        Network[] networksForTest = {s2bl1};

        NcmSimulator ncmSim = new NcmSimulator();

        for (Network net: networksForTest) {
            _prior_bl_mean = 1.0;
            _prior_bl_max = Double.POSITIVE_INFINITY;
            List<Tree> simulatedTrees = ncmSim.simulateTrees(net, num_settings, trees_per_setting);
//            System.out.print(simulatedTrees.toString());
            Map<Tree, Double> frequencies = ncmSim.getTopologyFrequencies(simulatedTrees);
//            System.out.println(frequencies.toString());

            gtpi.setBranchLengthExponentialPrior(_prior_bl_mean, _prior_bl_max);

            List<Tree> simulatedTopologies = new ArrayList<>();
            for (Tree t: frequencies.keySet()) {
                simulatedTopologies.add(t);
            }

            double mean_error = 0.0;
            double reference_error = 0.0;

            List<Double> probabilities = gtpi.calculateGTDistribution(net, simulatedTopologies, null, false);
            List<Double> pdf = normalizeDistribution(probabilities);
            List<Double> cdf = pdfToCdf(pdf);

            Map<Tree, Double> reference = new HashMap<>();
            int index;
            for (Tree top: simulatedTopologies) {
                reference.put(top, 0.0);
            }
            for (Tree st: simulatedTrees) {
                index = randomChoice(cdf);
                reference.put(simulatedTopologies.get(index), 1.0 / simulatedTrees.size() + reference.get(simulatedTopologies.get(index)));
            }

//            System.out.println(probabilities);
            int i = 0;
            for (Tree tree: simulatedTopologies) {
                System.out.println(tree.toString() + " freq: " + frequencies.get(tree) + ", prob: " + probabilities.get(i) + ", ref: " + reference.get(tree));
                double ratio = reference.get(tree) / probabilities.get(i);
                int num_samples = 0;

                if (ratio < 1) {
                    ratio = 1.0 / ratio;
                }
                if (ratio < Double.POSITIVE_INFINITY) {
                    reference_error += ratio;
                    System.out.print("reference error got " + ratio);
                } else {
                    reference_error += 1.0;
                    System.out.print("tree " + tree.toString() + " did not appear in reference sample");

                }

                ratio = frequencies.get(tree) / probabilities.get(i);
                if (ratio < 1) {
                    ratio = 1.0 / ratio;
                }
                mean_error += ratio;
                System.out.println(", freq error got " + ratio);
                i ++;

            }
            mean_error /= simulatedTopologies.size();
            reference_error /= simulatedTopologies.size();
            System.out.println("\nmean error: " + mean_error);
            System.out.println("reference error: " + reference_error);

        }
    }

    private void testBasicSampling() {
        List<Double> testPdf = new ArrayList<>();
        testPdf.add(1.);
        testPdf.add(2.);
        testPdf.add(3.);
        testPdf.add(7.);
        testPdf = normalizeDistribution(testPdf);
        System.out.println(testPdf.toString());
        List<Double> cdf = pdfToCdf(testPdf);
        System.out.println(cdf.toString());
        int i;
        for (i = 0; i < 20; i++) {
            System.out.println(i / 20.0 + ", " + binarySearch(cdf, i / 20.0));
        }
    }

    private List<Double> pdfToCdf(List<Double> pdf) {
        List<Double> cdf = new ArrayList<>(pdf.size());
        cdf.add(pdf.get(0));
        int i;
        for (i = 1; i < pdf.size(); i++) {
            cdf.add(cdf.get(i - 1) + pdf.get(i));
        }
        return cdf;
    }

    Random rand = new Random(314159);

    private int randomChoice(List<Double> cdf) {
        double r = rand.nextDouble();
        int choice = binarySearch(cdf, r);
        return choice;
    }

    private int binarySearch(List<Double> sorted, double value) {
        List<Double> sortedCopy = new ArrayList<>(sorted);
        while (sortedCopy.size() > 1) {
            int half = (sortedCopy.size() - 1) / 2;
            if (value > sortedCopy.get(half)) {
                sortedCopy = sortedCopy.subList(half + 1, sortedCopy.size());
            } else {
                sortedCopy = sortedCopy.subList(0, half + 1);
            }
        }
        return sorted.indexOf(sortedCopy.get(0));
    }

    private List<Double> normalizeDistribution(List<Double> dist) {
        List<Double> normalized = new ArrayList<>();
        double total = 0;
        for (Double d : dist) {
            total += d;
        }
        if (total <= 0) {
            System.err.print("Tried to normalize a distribution with sum <= 0");
            return null; // todo throw an exception
        }
        for (Double d : dist) {
            normalized.add(d / total);
        }
        return normalized;
    }
}
