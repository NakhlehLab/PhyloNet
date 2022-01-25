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
        String richNewick1 = "(E:10.0,((((B:2.0,(C:1.0)I7#H2:1.0::0.5)I6:2.0)I4#H1:2.0::0.5,(I7#H2:2.0::0.5,D:3.0)I5:3.0)I3:2.0,(A:5.0,I4#H1:1.0::0.5)I2:3.0)I1:2.0)I0;";

        scalingExperiment(richNewick1, 10, 1.5, 0.1);
        System.out.println();
        uniformScalingExperiment(richNewick1, 10, 1.2);
        System.out.println();
        reticulationExperiment(richNewick1, 0);
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
        Network<BniNetwork> newNetwork = originalNetwork.clone();

        System.out.println("Reticulation experiment (startReti=" + originalNetwork.getReticulationCount() +", minReti=" + minRet + ")");
        System.out.println("Reticulation Offset,Dissimilarity");


        for (int i = originalNetwork.getReticulationCount(); i >= minRet; i--) {
            TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity(originalNetwork, newNetwork);
            double dissimilarity = treeBasedDissimilarity.computeRootedBranchScore();

            System.out.println((i - originalNetwork.getReticulationCount()) + "," + dissimilarity);

            if (newNetwork.getReticulationCount() == 0)
                break;
            newNetwork.getNetworkNodes();
            ReticulationEdgeDeletion ed = new ReticulationEdgeDeletion();
            ed.setParameters(newNetwork, getReticulationEdges(newNetwork).get(0), null, null);
            ed.performOperation();
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
