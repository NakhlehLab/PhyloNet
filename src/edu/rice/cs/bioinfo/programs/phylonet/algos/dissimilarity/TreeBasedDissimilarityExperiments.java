package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.ReticulationEdgeDeletion;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.List;

public class TreeBasedDissimilarityExperiments {
    public static void main(String[] args) {
        if (args.length != 1)
            return;

        String richNewick = args[0];

        scalingExperiment(richNewick, 10, 1.5, 0.1);
        System.out.println();
        uniformScalingExperiment(richNewick, 10, 1.2);
        System.out.println();
        reticulationExperiment(richNewick, 0);
    }

    private static void uniformScalingExperiment(String richNewick, int numScales, double scaleFactor) {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);

        System.out.println("Uniform scaling experiment (numScales=" + numScales +", scaleFactor=" + scaleFactor + ")");
        System.out.println("Scale Factor,Dissimilarity");

        double curScale = 1.0;
        for (int i = 0; i <= numScales; i++) {
            Network<BniNetwork> newNetwork = originalNetwork.clone();
            Networks.scaleNetwork(newNetwork, curScale *= scaleFactor);

            TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity(originalNetwork, newNetwork);
            double dissimilarity = treeBasedDissimilarity.computeRootedBranchScore();

            System.out.println(curScale + "," + dissimilarity);
        }
    }

    private static void scalingExperiment(String richNewick, int numScales, double mean, double std) {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);
        Network<BniNetwork> currentNetwork = Networks.readNetwork(richNewick);

        NormalDistribution normalDistribution = new NormalDistribution(mean, std);

        System.out.println("Scaling experiment (numScales=" + numScales +", mean=" + mean + ", std=" + std + ")");
        System.out.println("# Iterations,Dissimilarity");

        for (int i = 0; i <= numScales; i++) {
            TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity(originalNetwork, currentNetwork);
            double dissimilarity = treeBasedDissimilarity.computeRootedBranchScore();

            System.out.println(i + "," + dissimilarity);

            for(NetNode<BniNetwork> node : Networks.postTraversal(currentNetwork)) {
                for(NetNode<BniNetwork> child : node.getChildren()) {
                    double scale = normalDistribution.sample();
                    child.setParentDistance(node, child.getParentDistance(node) * scale);
                }
            }
        }
    }

    private static void reticulationExperiment(String richNewick, int minRet) {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);
        Network<BniNetwork> currentNetwork = Networks.readNetwork(richNewick);

        System.out.println("Reticulation experiment (startReti=" + originalNetwork.getReticulationCount() +", minReti=" + minRet + ")");
        System.out.println("Reticulations Deleted,Dissimilarity");


        for (int i = originalNetwork.getReticulationCount(); i >= minRet; i--) {
            TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity(originalNetwork, currentNetwork);
            double dissimilarity = treeBasedDissimilarity.computeRootedBranchScore();

            System.out.println(originalNetwork.getReticulationCount() - i + "," + dissimilarity);

            if (currentNetwork.getReticulationCount() == 0)
                break;

            int j = 0;
            while (true) {
                ReticulationEdgeDeletion ed = new ReticulationEdgeDeletion();
                if (getReticulationEdges(currentNetwork).size() <= j)
                    throw new RuntimeException("Failed to perform reticulation deletion");

                ed.setParameters(currentNetwork, getReticulationEdges(currentNetwork).get(j), null, null);
                if (ed.performOperation())
                    break;

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
