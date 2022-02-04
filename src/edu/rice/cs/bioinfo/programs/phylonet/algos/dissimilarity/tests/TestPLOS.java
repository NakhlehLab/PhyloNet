package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.tests;

import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.R;
import edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.RootedNetworkBranchScore;
import edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.WeightedAveragePathDistance;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class TestPLOS {
    static List<Double> lambdas = new ArrayList<>();

    static {
        lambdas.add(0.0);
        lambdas.add(1.0);
        lambdas.add(3.0);
        lambdas.add(2.0);
        lambdas.add(5.0);
        lambdas.add(4.0);
        lambdas.add(2.0);
        lambdas.add(7.0);
        lambdas.add(4.0);
        lambdas.add(6.0);
        lambdas.add(4.0);
        lambdas.add(8.0);
        lambdas.add(2.0);
    }

    private static Network<BniNetwork> getNetwork1() {
        Network<BniNetwork> plosNet = Networks.readNetwork("Root;");

        NetNode<BniNetwork> root = plosNet.getRoot();
        NetNode<BniNetwork> a = new BniNetNode<>();
        a.setName("a");
        NetNode<BniNetwork> b = new BniNetNode<>();
        b.setName("b");
        NetNode<BniNetwork> c = new BniNetNode<>();
        c.setName("c");
        NetNode<BniNetwork> d = new BniNetNode<>();
        d.setName("d");

        NetNode<BniNetwork> I1 = new BniNetNode<>();
        NetNode<BniNetwork> I2 = new BniNetNode<>();
        NetNode<BniNetwork> I3 = new BniNetNode<>();
        NetNode<BniNetwork> I4 = new BniNetNode<>();
        NetNode<BniNetwork> I5 = new BniNetNode<>();
        NetNode<BniNetwork> I6 = new BniNetNode<>();

        root.adoptChild(I1, lambdas.get(1));
        root.adoptChild(I2, lambdas.get(2));

        I1.adoptChild(a, lambdas.get(12));
        I1.adoptChild(I6, lambdas.get(6));

        I2.adoptChild(I3, lambdas.get(3));
        I2.adoptChild(I4, lambdas.get(10));

        I6.adoptChild(b, lambdas.get(7));

        I3.adoptChild(d, lambdas.get(4));
        I3.adoptChild(I5, lambdas.get(5));

        I4.adoptChild(c, lambdas.get(11));
        I4.adoptChild(I5, lambdas.get(9));

        I5.adoptChild(I6, lambdas.get(8));

        return plosNet;
    }

    private static Network<BniNetwork> getNetwork2() {
        Random random = new Random();
        double y = random.nextDouble() * lambdas.get(7);
        double x = -y + (Math.min(lambdas.get(6), lambdas.get(5) + lambdas.get(8)) + y) * random.nextDouble();

        Network<BniNetwork> plosNet = Networks.readNetwork("Root;");

        NetNode<BniNetwork> root = plosNet.getRoot();
        NetNode<BniNetwork> a = new BniNetNode<>();
        a.setName("a");
        NetNode<BniNetwork> b = new BniNetNode<>();
        b.setName("b");
        NetNode<BniNetwork> c = new BniNetNode<>();
        c.setName("c");
        NetNode<BniNetwork> d = new BniNetNode<>();
        d.setName("d");

        NetNode<BniNetwork> I1 = new BniNetNode<>();
        NetNode<BniNetwork> I2 = new BniNetNode<>();
        NetNode<BniNetwork> I3 = new BniNetNode<>();
        NetNode<BniNetwork> I4 = new BniNetNode<>();
        NetNode<BniNetwork> I5 = new BniNetNode<>();
        NetNode<BniNetwork> I6 = new BniNetNode<>();

        root.adoptChild(I1, lambdas.get(1));
        root.adoptChild(I2, lambdas.get(2));

        I1.adoptChild(a, lambdas.get(12));
        I1.adoptChild(I5, lambdas.get(6) - x);

        I2.adoptChild(I3, lambdas.get(3));
        I2.adoptChild(I4, lambdas.get(10));

        I5.adoptChild(I6, x + y);

        I3.adoptChild(d, lambdas.get(4));
        I3.adoptChild(I5, lambdas.get(5) + lambdas.get(8) - x);

        I4.adoptChild(c, lambdas.get(11));
        I4.adoptChild(I6, lambdas.get(8) + lambdas.get(9) + y);

        I6.adoptChild(b, lambdas.get(7) - y);

        return plosNet;
    }

    public static void main(String[] args) {
        Network<BniNetwork> net1 = getNetwork1();
        Network<BniNetwork> net2 = getNetwork2();

        System.out.println(net1);
        for (NetworkTree<BniNetwork> tree : Networks.getTrees(net1)) {
            System.out.println(tree.makeTree());
        }

        System.out.println();

        System.out.println(net2);
        for (NetworkTree<BniNetwork> tree : Networks.getTrees(net2)) {
            System.out.println(tree.makeTree());
        }

        System.out.println();

        RootedNetworkBranchScore RNBS = new RootedNetworkBranchScore(net1, net2);
        System.out.println("rNBS: " + RNBS.compute());
        System.out.println("NormRNBS: " + RNBS.computeNormalized());

        WeightedAveragePathDistance WAPD = new WeightedAveragePathDistance(net1, net2);
        System.out.println("WAPD: " + WAPD.compute());
        System.out.println("NormWAPD: " + WAPD.computeNormalized());

    }
}
