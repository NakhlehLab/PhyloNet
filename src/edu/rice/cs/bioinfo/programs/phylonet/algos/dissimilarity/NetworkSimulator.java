package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class NetworkSimulator <T> {
    private final Random random;
    private final Network<T> originalNetwork;

    public NetworkSimulator(int seed) {
        this.random = new Random(seed);
        this.originalNetwork = getNewNetwork();

    }

    public NetworkSimulator() {
        this.random = new Random();
        this.originalNetwork  = getNewNetwork();

        System.out.println(this.originalNetwork.toString());
    }

    private Network<T> generateNewRootedNetwork(List<String> leafSet) {
        Network<T> network = Networks.readNetwork("Root;");

        return network;
    }

    private Network<T> getNewNetwork() {
        List<Double> blPool = new ArrayList<>();

        for (int i = 0; i < 12; i++)
            blPool.add(0.7 + this.random.nextDouble() * (1.3 - 0.7));

        double blTotal = 8.0;
        double blI1 = 2.0;
        double blJ1 = 2.0;
        double blX = blPool.get(0);
        double blY = blPool.get(1);
        double blS = blPool.get(2);
        double blT = blPool.get(3);
        double blV = blPool.get(4);
        double blW = blPool.get(5);
        double blO = blPool.get(6);
        double blP = blPool.get(7);
        double blQ = blPool.get(8);
        double blR = blPool.get(9);
        double prob1 = 0.5;
        double prob2 = 1 - prob1;
        double prob3 = 0.5;
        double prob4 = 1 - prob3;


        Network<T> trueNetwork = Networks.readNetwork("Root;");

        NetNode<T> nodeX = new BniNetNode<>();
        trueNetwork.getRoot().adoptChild(nodeX, blX);

        NetNode<T> nodeY = new BniNetNode<>();
        trueNetwork.getRoot().adoptChild(nodeY, blY);

        double sumOP = (blS + blT + blO + blP) / 2;
        double sumO = blS + blO;
        double sumP = blT + blP;
        blS = blS / sumO * sumOP;
        blO = blO / sumO * sumOP;
        blT = blT / sumP * sumOP;
        blP = blP / sumP * sumOP;

        NetNode<T> nodeS = new BniNetNode<>();
        nodeX.adoptChild(nodeS, blS);
        NetNode<T> nodeK = new BniNetNode<>();
        nodeK.setName("K");
        nodeS.adoptChild(nodeK, blTotal - blS - blX);
        NetNode<T> nodeO = new BniNetNode<>();
        nodeS.adoptChild(nodeO, blO);


        NetNode<T> nodeT = new BniNetNode<>();
        nodeX.adoptChild(nodeT, blT);
        NetNode<T> nodeL = new BniNetNode<>();
        nodeL.setName("L");
        nodeT.adoptChild(nodeL, blTotal - blT - blX);
        NetNode<T> nodeP = new BniNetNode<>();
        nodeT.adoptChild(nodeP, blP);

        NetNode<T> nodeI = new BniNetNode<>();
        nodeO.adoptChild(nodeI, 0.0);
        nodeP.adoptChild(nodeI, 0.0);
        nodeI.setParentProbability(nodeO, prob1);
        nodeI.setParentProbability(nodeP, prob2);
        NetNode<T> nodeI1 = new BniNetNode<>();
        nodeI.adoptChild(nodeI1, blI1);

        NetNode<T> nodeA = new BniNetNode<>();
        nodeA.setName("A");
        NetNode<T> nodeB = new BniNetNode<>();
        nodeB.setName("B");
        NetNode<T> nodeC = new BniNetNode<>();
        nodeC.setName("C");
        NetNode<T> nodeD = new BniNetNode<>();
        nodeD.setName("D");

        nodeO.adoptChild(nodeA, blTotal - sumOP - blX);
        nodeI1.adoptChild(nodeB, blTotal - sumOP - blX - blI1);
        nodeI1.adoptChild(nodeC, blTotal - sumOP - blX - blI1);
        nodeP.adoptChild(nodeD, blTotal - sumOP - blX);




        double sumQR = (blV + blQ + blW + blR) / 2;
        double sumQ = blV + blQ;
        double sumR = blW + blR;
        blV = blV / sumQ * sumQR;
        blQ = blQ / sumQ * sumQR;
        blW = blW / sumR * sumQR;
        blR = blR / sumR * sumQR;

        NetNode<T> nodeV = new BniNetNode<>();
        nodeY.adoptChild(nodeV, blV);
        NetNode<T> nodeM = new BniNetNode<>();
        nodeM.setName("M");
        nodeV.adoptChild(nodeM, blTotal - blV - blY);
        NetNode<T> nodeQ = new BniNetNode<>();
        nodeV.adoptChild(nodeQ, blQ);


        NetNode<T> nodeW = new BniNetNode<>();
        nodeY.adoptChild(nodeW, blW);
        NetNode<T> nodeN = new BniNetNode<>();
        nodeN.setName("N");
        nodeW.adoptChild(nodeN, blTotal - blW - blY);
        NetNode<T> nodeR = new BniNetNode<>();
        nodeW.adoptChild(nodeR, blR);

        NetNode<T> nodeJ = new BniNetNode<>();
        nodeQ.adoptChild(nodeJ, 0.0);
        nodeR.adoptChild(nodeJ, 0.0);
        nodeJ.setParentProbability(nodeQ, prob3);
        nodeJ.setParentProbability(nodeR, prob4);
        NetNode<T> nodeJ1 = new BniNetNode<>();
        nodeJ.adoptChild(nodeJ1, blJ1);

        NetNode<T> nodeE = new BniNetNode<>();
        nodeE.setName("E");
        NetNode<T> nodeF = new BniNetNode<>();
        nodeF.setName("F");
        NetNode<T> nodeG = new BniNetNode<>();
        nodeG.setName("G");
        NetNode<T> nodeH = new BniNetNode<>();
        nodeH.setName("H");

        nodeQ.adoptChild(nodeE, blTotal - sumQR - blY);
        nodeJ1.adoptChild(nodeF, blTotal - sumQR - blY - blJ1);
        nodeJ1.adoptChild(nodeG, blTotal - sumQR - blY - blJ1);
        nodeR.adoptChild(nodeH, blTotal - sumQR - blY);



        return trueNetwork;
    }
}
