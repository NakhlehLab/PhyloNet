package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

public class MatrixBasedDissimilarityExperiments {
    public static void main(String[] args) {
        String richNewick1 = "(E:10.0,((((B:2.0,(C:1.0)I7#H2:1.0::0.5)I6:2.0)I4#H1:2.0::0.5,(I7#H2:2.0::0.5,D:3.0)I5:3.0)I3:2.0,(A:5.0,I4#H1:1.0::0.5)I2:3.0)I1:2.0)I0;";

        uniformScalingExperiment(richNewick1, 10, 1.2);
    }

    private static void uniformScalingExperiment(String richNewick, int numScales, double scaleFactor) {
        Network<BniNetwork> originalNetwork = Networks.readNetwork(richNewick);

        double curScale = 1.0;
        for (int i = 0; i <= numScales; i++) {
            Network<BniNetwork> newNetwork = originalNetwork.clone();
            Networks.scaleNetwork(newNetwork, curScale *= scaleFactor);

            MatrixBasedDissimilarity<BniNetwork> matrixBasedDissimilarity = new MatrixBasedDissimilarity(originalNetwork, newNetwork);
            double dissimilarity = matrixBasedDissimilarity.computeMeanDistance();

//            System.out.println("Scale Factor: " + curScale + " - Dissimilarity: " + dissimilarity);
            System.out.println(curScale + "," + dissimilarity);
        }
    }
}
