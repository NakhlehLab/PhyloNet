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

public class RuntimeExperiments {
    public static Map<String, BiFunction<Network<BniNetwork>, Network<BniNetwork>, Double>> metrics = new HashMap<>();
    static {
        metrics.put("WAPD", RuntimeExperiments::getWAPD);
        metrics.put("normWAPD", RuntimeExperiments::getNormWAPD);
        metrics.put("APD", RuntimeExperiments::getAPD);
        metrics.put("normAPD", RuntimeExperiments::getNormAPD);
        metrics.put("rNBS", RuntimeExperiments::getRNBS);
        metrics.put("NormRNBS", RuntimeExperiments::getNormRNBS);
        metrics.put("Nakhleh", RuntimeExperiments::getNakhleh);
    }

    private static String metric = "WAPD";
    private static int maxReti = 10;
    private static Network<BniNetwork> baseNet = Networks.readNetwork("(((((E:0.006688512350679896)#H1[&gamma=0.49898092591168197]:0.01847451738763451,(#H1:0.011148402346383788,(D:0.010623666508974701,((G:0.002022549370932501,H:0.002022549370932501)S9:0.0024288541668802035,F:0.004451403537812705)S8:0.006172262971161997)S7:0.007213248188088982)S6:0.007326115041250723)S5:0.006963239804625772,(C:0.017973898313667334)#H2[&gamma=0.5864355936302258]:0.014152371229272844)S2:0.009136331956569717,((B:0.02189207431838481,#H2:0.003918176004717475)S4:0.014810291498616228,A:0.03670236581700104)S3:0.004560235682508858)S1:0.05873739850049011);");

    public static void main(String[] args) {
        if (args.length != 3) {
            System.err.println("Command-line arguments:\n" +
                    "\t[WAPD/normWAPD/rNBS/NormRNBS/Nakhleh] [.trees file] [results directory]");
            System.exit(-1);
        }

        metric = args[0];

        if (!metrics.containsKey(metric)) {
            System.err.println("Can use the following measures: WAPD/normWAPD/rNBS/NormRNBS/Nakhleh");
            System.exit(-1);
        }

        List<Network<BniNetwork>> networks = new ArrayList<>();

        try {
            Scanner fileReader = new Scanner(new File(args[1]));

            while (fileReader.hasNextLine()) {
                String line = fileReader.nextLine();

                if (line.startsWith("tree ")) {
                    networks.add(Networks.readNetwork(line.substring(line.indexOf("=") + 1)));
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

        try {
            BufferedWriter edgeCountWriter = new BufferedWriter(new FileWriter(Paths.get(args[2], "edgecount_runtime.txt").toFile()));
            BufferedWriter retCountWriter = new BufferedWriter(new FileWriter(Paths.get(args[2], "retcount_runtime.txt").toFile()));

            edgeCountWriter.write("Runtime vs. number of edges experiment\n");
            edgeCountWriter.write("Edge Count,Runtime (ms)\n");

            retCountWriter.write("Runtime vs. number of reticulations experiment\n");
            retCountWriter.write("Reticulation Count,Runtime (ms)\n");

            for (int i = 0; i < networks.size(); i++) {
                if (networks.get(i).getReticulationCount() > maxReti) {
                    System.out.println("Skipping " + i);
                    continue;
                } else {
                    System.out.println("Running " + i);
                }

                long startTime = System.nanoTime();
                double dissimilarity = metrics.get(metric).apply(baseNet, networks.get(i));
                long endTime = System.nanoTime();
                double runtime = (endTime - startTime) / 1000000.0;

                edgeCountWriter.write(networks.get(i).getEdgeCount() + "," + runtime + "\n");
                retCountWriter.write(networks.get(i).getReticulationCount() + "," + runtime + "\n");
            }

            edgeCountWriter.close();
            retCountWriter.close();

        } catch (IOException e) {
            System.err.println("Could not write experiment results.");
            e.printStackTrace();
            System.exit(-1);
        }
    }

    private static double getWAPD(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
        return WeightedAveragePathDistance.compute(originalNetwork, currentNetwork, true);
    }

    private static double getNormWAPD(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
        return WeightedAveragePathDistance.computeNormalized(originalNetwork, currentNetwork, true);
    }

    private static double getAPD(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
        return AveragePathDistance.compute(originalNetwork, currentNetwork);
    }

    private static double getNormAPD(Network<BniNetwork> originalNetwork, Network<BniNetwork> currentNetwork) {
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
