package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.experiments;

import edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.AveragePathDistance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.RootedNetworkBranchScore;
import edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.WeightedAveragePathDistance;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.BiFunction;

public class McmcseqPairwiseExperiments {

    public static Map<String, BiFunction<Network<BniNetwork>, Network<BniNetwork>, Double>> metrics = new HashMap<>();
    public static Map<String, BiFunction<double[][], double[][], Double>> metrics2 = new HashMap<>();
    static {
        metrics2.put("WAPD", McmcseqPairwiseExperiments::getWAPD);
        metrics2.put("normWAPD", McmcseqPairwiseExperiments::getNormWAPD);
        metrics2.put("APD", McmcseqPairwiseExperiments::getAPD);
        metrics2.put("normAPD", McmcseqPairwiseExperiments::getNormAPD);
        metrics.put("rNBS", McmcseqPairwiseExperiments::getRNBS);
        metrics.put("NormRNBS", McmcseqPairwiseExperiments::getNormRNBS);
        metrics.put("Nakhleh", McmcseqPairwiseExperiments::getNakhleh);
    }

    private static String metric = "WAPD";

    public static void main(String[] args) {
        if (args.length != 4) {
            System.err.println("Command-line arguments:\n" +
                    "\t[WAPD/normWAPD/rNBS/NormRNBS/Nakhleh] [.out file] [results directory] [Nth Network]");
            System.exit(-1);
        }

        metric = args[0];

        if (!metrics.containsKey(metric) && !metrics2.containsKey(metric)) {
            System.err.println("Can use the following measures: WAPD/normWAPD/rNBS/NormRNBS/Nakhleh");
            System.exit(-1);
        }

        List<Network<BniNetwork>> networks = new ArrayList<>();
        List<double[][]> networksMatrices = new ArrayList<>();

        try {
            Scanner fileReader = new Scanner(new File(args[1]));

            while (fileReader.hasNextLine()) {
                String line = fileReader.nextLine();

                if (line.startsWith("[")) {
                    String richNewick = line.substring(line.indexOf("]") + 1);
                    networks.add(Networks.readNetwork(richNewick));
                    networksMatrices.add(WeightedAveragePathDistance.computeMatrix(Networks.readNetwork(richNewick), true));
                }
            }
            fileReader.close();

        } catch (FileNotFoundException e) {
            System.err.println("Could not find .out file to read MCMC_SEQ posterior samples from.");
            System.exit(-1);
        }

        File saveDir = new File(args[2]);

        if (!saveDir.exists() || !saveDir.isDirectory()) {
            System.err.println("The results directory to write experiment results to does not exist.");
            System.exit(-1);
        }

        int Nth = Integer.parseInt(args[3]);

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(Paths.get(args[2], "mcmcseq_pairwise.txt").toFile()));

            writer.write("MCMC_SEQ experiment\n");
            writer.write("X,Y,Z\n");

            for (int i = 0; i < networks.size(); i++) {
                if (i % Nth != 0)
                    continue;
                for (int j = 0; j <= i; j++) {
                    if (j % Nth != 0)
                        continue;

                    double dissimilarity = Double.NaN;
                    if (metrics.containsKey(metric))
                        dissimilarity = metrics.get(metric).apply(networks.get(i), networks.get(j));
                    else if (metrics2.containsKey(metric))
                        dissimilarity = metrics2.get(metric).apply(networksMatrices.get(i), networksMatrices.get(j));

                    writer.write(i + "," + j + "," + dissimilarity + "\n");
                }
                System.out.println(i + "/" + (networks.size() - 1));
            }

            writer.close();

        } catch (IOException e) {
            System.err.println("Could not write experiment results.");
            e.printStackTrace();
            System.exit(-1);
        }
    }

    private static double getWAPD(double[][] originalNetwork, double[][] currentNetwork) {
        return WeightedAveragePathDistance.compute(originalNetwork, currentNetwork);
    }

    private static double getNormWAPD(double[][] originalNetwork, double[][] currentNetwork) {
        return WeightedAveragePathDistance.computeNormalized(originalNetwork, currentNetwork);
    }

    private static double getAPD(double[][] originalNetwork, double[][] currentNetwork) {
        return AveragePathDistance.compute(originalNetwork, currentNetwork);
    }

    private static double getNormAPD(double[][] originalNetwork, double[][] currentNetwork) {
        return AveragePathDistance.computeNormalized(originalNetwork, currentNetwork);
    }

    private static double getRNBS(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
        return RootedNetworkBranchScore.compute(originalNetwork, currentNetwork);
    }

    private static double getNormRNBS(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
        return RootedNetworkBranchScore.computeNormalized(originalNetwork, currentNetwork);
    }

    private static double getNakhleh(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
        return Networks.computeDistanceBetweenTwoNetworks(originalNetwork, currentNetwork);
    }
}
