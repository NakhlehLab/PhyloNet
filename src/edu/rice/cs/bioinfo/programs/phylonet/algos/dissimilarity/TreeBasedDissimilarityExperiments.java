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
        String richNewick2 = "((((((F:0.013882235548125541)#H3[&gamma=0.013473703118829228]:0.0034580340447386143,E:0.017340269592864156)S10:0.023538514552048823)#H6[&gamma=0.04622625080145859]:0.024057416697191088,(#H6:0.009797503385167383)#H7[&gamma=0.8327360978123327]:0.014259913312023705)S2:0.012551291969394363,((#H7:0.0071200282766816225,((((H:0.0019152517083881808)#H1[&gamma=0.228212416050242]:0.016597481307807815)#H4[&gamma=0.32745258490120244]:0.012408862149908363,(#H3:0.011318078666590185,C:0.025200314214715726)S13:0.005721280951388633)S7:0.022053398150367293,((D:0.022562929979630955)#H5[&gamma=0.936402574589784]:0.018528227174561165,((#H5:0.00703283054045933,#H4:0.01108302750389429)S11:0.007006486390104139,((#H1:0.007961293406641073,((G:0.003864190533995076)#H2[&gamma=0.8111074321492627]:0.003612132656480802,#H2:0.003612132656480802)S15:0.002400221924553376)S14:0.015329554116387891,B:0.025206099231417145)S12:0.01139614767877728)S9:0.004488910243997696)S8:0.011883836162279532)S6:0.004821322490290332)S5:0.015214760220033062,((A:0.05352310311604339)#H8[&gamma=0.13029728722782097]:7.49918073490792E-4,#H8:7.49918073490792E-4)S4:0.018738054837260862)S3:0.0044764167847033826)S1:0.022512507188501577);";
        scalingExperiment(richNewick2, 10, 1.5, 0.1);
        System.out.println();
        uniformScalingExperiment(richNewick2, 10, 1.2);
        System.out.println();
        reticulationExperiment(richNewick2, 0);
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
        System.out.println("Reticulation Offset,Dissimilarity");


        for (int i = originalNetwork.getReticulationCount(); i >= minRet; i--) {
            TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity(originalNetwork, currentNetwork);
            double dissimilarity = treeBasedDissimilarity.computeRootedBranchScore();

            System.out.println((i - originalNetwork.getReticulationCount()) + "," + dissimilarity);

            if (currentNetwork.getReticulationCount() == 0)
                break;

            ReticulationEdgeDeletion ed = new ReticulationEdgeDeletion();
            ed.setParameters(currentNetwork, getReticulationEdges(currentNetwork).get(0), null, null);
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
