package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

public class DissimilarityExperiments {
    public static void main(String[] args) {
        String richNewick1 = "(E:10.0,((((B:2.0,(C:1.0)I7#H2:1.0::0.5)I6:2.0)I4#H1:2.0::0.5,(I7#H2:2.0::0.5,D:3.0)I5:3.0)I3:2.0,(A:5.0,I4#H1:1.0::0.5)I2:3.0)I1:2.0)I0;";

        reticulationExperiment(richNewick1, 0, 5);
    }

    private static void reticulationExperiment(String richNewick, int minRet, int maxRet) {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);

        for (int i = minRet; i <= maxRet; i++) {
            Network<BniNetwork> newNetwork = originalNetwork.clone();
            Networks.addRandomReticulationEdge(newNetwork, i);
            Networks.autoLabelNodes(newNetwork);

            MatrixBasedDissimilarity<BniNetwork> matrixBasedDissimilarity = new MatrixBasedDissimilarity(originalNetwork, newNetwork);

            double dissimilarity = matrixBasedDissimilarity.computeMeanDistance();

            System.out.println("# Added reticulations: " + i + " - Dissimilarity: " + dissimilarity);
        }
    }
}
