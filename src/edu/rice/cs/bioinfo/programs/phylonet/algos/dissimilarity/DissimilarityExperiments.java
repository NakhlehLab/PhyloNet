package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

public class DissimilarityExperiments {
    public static void main(String[] args) {
        String richNewick1 = "(E:10.0,((((B:2.0,(C:1.0)I7#H2:1.0::0.5)I6:2.0)I4#H1:2.0::0.5,(I7#H2:2.0::0.5,D:3.0)I5:3.0)I3:2.0,(A:5.0,I4#H1:1.0::0.5)I2:3.0)I1:2.0)I0;";

        uniformScalingExperiment(richNewick1, 5, 1.2);
    }

    private static void uniformScalingExperiment(String richNewick, int numScales, double scaleFactor) {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);

        double curScale = 1.0;
        for (int i = 0; i <= numScales; i++) {
            Network<BniNetwork> newNetwork = originalNetwork.clone();
            Networks.scaleNetwork(newNetwork, curScale *= scaleFactor);

            TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity(originalNetwork, newNetwork);
            double dissimilarity = treeBasedDissimilarity.computeRootedBranchScore();

            System.out.println("Scale Factor: " + curScale + " - Dissimilarity: " + dissimilarity);
        }
    }

    private static void reticulationExperiment(String richNewick, int minRet, int maxRet) {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);

        for (int i = minRet; i <= maxRet; i++) {
            Network<BniNetwork> newNetwork = originalNetwork.clone();
            Networks.addRandomReticulationEdge(newNetwork, i);
            Networks.autoLabelNodes(newNetwork);


            TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity(originalNetwork, newNetwork);
            double dissimilarity = treeBasedDissimilarity.computeRootedBranchScore();

            System.out.println("# Added reticulations: " + i + " - Dissimilarity: " + dissimilarity);
        }
    }
}
