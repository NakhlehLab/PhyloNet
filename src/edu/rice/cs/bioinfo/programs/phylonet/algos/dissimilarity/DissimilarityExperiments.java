package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.ReticulationEdgeDeletion;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.*;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.BiFunction;

public class DissimilarityExperiments {

    public static Map<String, BiFunction<Network<BniNetwork>, Network<BniNetwork>, Double>> metrics = new HashMap<>();
    static {
        metrics.put("WAPD", DissimilarityExperiments::getWAPD);
        metrics.put("normWAPD", DissimilarityExperiments::getNormWAPD);
        metrics.put("rNBS", DissimilarityExperiments::getRNBS);
        metrics.put("NormRNBS", DissimilarityExperiments::getNormRNBS);
    }

    private static String metric = "WAPD";
    private static int maxReti = 10;

    public static void main(String[] args) {
        if (args.length != 3) {
            System.err.println("Command-line arguments:\n" +
                    "\t[WAPD/normWAPD/rNBS/NormRNBS] [.trees file] [results directory]");
            System.exit(-1);
        }

        metric = args[0];

        if (!metrics.containsKey(metric)) {
            System.err.println("Can use the following measures: WAPD/normWAPD/rNBS/NormRNBS");
            System.exit(-1);
        }

        List<String> richNewicks = new ArrayList<>();

        try {
            Scanner fileReader = new Scanner(new File(args[1]));

            while (fileReader.hasNextLine()) {
                String line = fileReader.nextLine();

                if (line.startsWith("tree ")) {
                    richNewicks.add(line.substring(line.indexOf("=") + 1));
                }
            }
            fileReader.close();

        } catch (FileNotFoundException e) {
            System.err.println("Could not find .trees file to read networks from.");
            System.exit(-1);
        }

        File saveDir = new File(args[2]);

        if (!saveDir.exists() || !saveDir.isDirectory()) {
            System.err.println("The results directory to write experiment results to does not exist.");
            System.exit(-1);
        }

        try {
            int i = 0;
            for (String richNewick : richNewicks) {

                if (Networks.readNetwork(richNewick).getReticulationCount() > maxReti) {
                    System.out.println("Skipping " + i);
                    continue;
                } else {
                    System.out.println("Running " + i);
                }

                BufferedWriter writer = new BufferedWriter(new FileWriter(Paths.get(args[2], "net" + i + ".txt").toFile()));

                scalingExperiment(richNewick, writer, 10, 1.5, 0.1);
                uniformScalingExperiment(richNewick, writer, 10, 1.2);
                reticulationExperiment(richNewick, writer, 0);

                i += 1;

                writer.close();
            }
        } catch (IOException e) {
            System.err.println("Could not write experiment results.");
            e.printStackTrace();
            System.exit(-1);
        }
    }

    private static void uniformScalingExperiment(String richNewick, BufferedWriter writer, int numScales, double scaleFactor) throws IOException {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);

        writer.write("Uniform scaling experiment (numScales=" + numScales +", scaleFactor=" + scaleFactor + ")\n");
        writer.write("Scale Factor,Dissimilarity\n");

        double curScale = 1.0;
        for (int i = 0; i <= numScales; i++) {
            Network<BniNetwork> newNetwork = originalNetwork.clone();
            Networks.scaleNetwork(newNetwork, curScale *= scaleFactor);

            double dissimilarity = metrics.get(metric).apply(originalNetwork, newNetwork);

            writer.write(curScale + "," + dissimilarity + "\n");
        }
    }

    private static void scalingExperiment(String richNewick, BufferedWriter writer, int numScales, double mean, double std) throws IOException {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);
        Network<BniNetwork> currentNetwork = Networks.readNetwork(richNewick);

        NormalDistribution normalDistribution = new NormalDistribution(mean, std);

        writer.write("Scaling experiment (numScales=" + numScales +", mean=" + mean + ", std=" + std + ")\n");
        writer.write("Number of Iterations,Dissimilarity\n");

        for (int i = 0; i <= numScales; i++) {
            double dissimilarity = metrics.get(metric).apply(originalNetwork, currentNetwork);

            writer.write(i + "," + dissimilarity + "\n");

            for(NetNode<BniNetwork> node : Networks.postTraversal(currentNetwork)) {
                for(NetNode<BniNetwork> child : node.getChildren()) {
                    double scale = normalDistribution.sample();
                    child.setParentDistance(node, child.getParentDistance(node) * scale);
                }
            }
        }
    }

    private static void reticulationExperiment(String richNewick, BufferedWriter writer, int minRet) throws IOException {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);
        Network<BniNetwork> currentNetwork = Networks.readNetwork(richNewick);

        writer.write("Reticulation experiment (startReti=" + originalNetwork.getReticulationCount() +", minReti=" + minRet + ")\n");
        writer.write("Number of Iterations,Dissimilarity\n");


        for (int i = originalNetwork.getReticulationCount(); i >= minRet; i--) {
            double dissimilarity = metrics.get(metric).apply(originalNetwork, currentNetwork);

            writer.write(originalNetwork.getReticulationCount() - i + "," + dissimilarity + "\n");

            if (currentNetwork.getReticulationCount() == 0)
                break;

            int j = 0;
            while (true) {
                ReticulationEdgeDeletion ed = new ReticulationEdgeDeletion();
                if (getReticulationEdges(currentNetwork).size() <= j) {
                    System.err.println("Failed to perform reticulation deletion");
                    break;
                }

                ed.setParameters(currentNetwork, getReticulationEdges(currentNetwork).get(j), null, null);
                try {
                    if (ed.performOperation())
                        break;
                } catch (RuntimeException e) {
                    System.err.println("Failed to perform reticulation deletion");
                    break;
                }

                j += 1;
            }
        }
    }

    private static List<Tuple<NetNode, NetNode>> getReticulationEdges(Network<BniNetwork> network) {
        List<Tuple<NetNode, NetNode>> list = new ArrayList<>();
        for (NetNode<BniNetwork> node : Networks.getInternalNodes(network)) {
            for (NetNode<BniNetwork> childNode : node.getChildren()) {
                if (node.isTreeNode() && childNode.isNetworkNode()) {
                    list.add(new Tuple<>(node , childNode));
                }
            }
        }

        return list;
    }

    private static double getWAPD(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
        WeightedAveragePathDistance<BniNetwork> WAPD = new WeightedAveragePathDistance(originalNetwork, currentNetwork);
        return WAPD.compute();
    }

    private static double getNormWAPD(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
        WeightedAveragePathDistance<BniNetwork> WAPD = new WeightedAveragePathDistance(originalNetwork, currentNetwork);
        return WAPD.computeNormalized();
    }

    private static double getRNBS(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
        RootedNetworkBranchScore<BniNetwork> RNBS = new RootedNetworkBranchScore(originalNetwork, currentNetwork);
        return RNBS.compute();
    }

    private static double getNormRNBS(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
        RootedNetworkBranchScore<BniNetwork> RNBS = new RootedNetworkBranchScore(originalNetwork, currentNetwork);
        return RNBS.computeNormalized();
    }
}
