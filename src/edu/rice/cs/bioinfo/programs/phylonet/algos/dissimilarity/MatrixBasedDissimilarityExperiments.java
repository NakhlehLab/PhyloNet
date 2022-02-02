package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.ReticulationEdgeDeletion;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class MatrixBasedDissimilarityExperiments {
    public static void main(String[] args) throws IOException {
        if (args.length != 2)
            return;

        List<String> richNewicks = new ArrayList<>();

        Scanner fileReader = new Scanner(new File(args[0]));
        while (fileReader.hasNextLine()) {
            String line = fileReader.nextLine();

            if (line.startsWith("tree ")) {
                richNewicks.add(line.substring(line.indexOf("=") + 1));
            }
        }
        fileReader.close();

        int i = 0;
        for (String richNewick : richNewicks) {
            System.out.println("Running " + i);
            BufferedWriter writer = new BufferedWriter(new FileWriter(args[1] + "/net" + i + ".txt"));

            scalingExperiment(richNewick, writer, 10, 1.5, 0.1);
            uniformScalingExperiment(richNewick, writer, 10, 1.2);
            reticulationExperiment(richNewick, writer, Math.max(Networks.readNetwork(richNewick).getReticulationCount() - 10, 0));

            i += 1;

            writer.close();
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

            MatrixBasedDissimilarity<BniNetwork> matrixBasedDissimilarity = new MatrixBasedDissimilarity(originalNetwork, newNetwork);
            double dissimilarity = matrixBasedDissimilarity.computeWeightedAveragePathDistance();

            writer.write(curScale + "," + dissimilarity + "\n");
        }
    }

    private static void scalingExperiment(String richNewick, BufferedWriter writer, int numScales, double mean, double std) throws IOException {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);
        Network<BniNetwork> currentNetwork = Networks.readNetwork(richNewick);

        NormalDistribution normalDistribution = new NormalDistribution(mean, std);

        writer.write("Scaling experiment (numScales=" + numScales +", mean=" + mean + ", std=" + std + ")\n");
        writer.write("# Iterations,Dissimilarity\n");

        for (int i = 0; i <= numScales; i++) {
            MatrixBasedDissimilarity<BniNetwork> matrixBasedDissimilarity = new MatrixBasedDissimilarity(originalNetwork, currentNetwork);
            double dissimilarity = matrixBasedDissimilarity.computeWeightedAveragePathDistance();

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
        writer.write("Reticulations Deleted,Dissimilarity\n");


        for (int i = originalNetwork.getReticulationCount(); i >= minRet; i--) {
            MatrixBasedDissimilarity<BniNetwork> matrixBasedDissimilarity = new MatrixBasedDissimilarity(originalNetwork, currentNetwork);
            double dissimilarity = matrixBasedDissimilarity.computeWeightedAveragePathDistance();

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
}
