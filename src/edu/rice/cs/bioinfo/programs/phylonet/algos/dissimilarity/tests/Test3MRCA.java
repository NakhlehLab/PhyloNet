package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.tests;

import edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.AveragePathDistance;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.Arrays;

public class Test3MRCA {
    private static Network<BniNetwork> getNetwork1() {
        Network<BniNetwork> testNet = Networks.readNetwork("Root;");

        NetNode<BniNetwork> root = testNet.getRoot();
        NetNode<BniNetwork> a = new BniNetNode<>();
        a.setName("a");
        NetNode<BniNetwork> b = new BniNetNode<>();
        b.setName("b");

        NetNode<BniNetwork> I1 = new BniNetNode<>();
        NetNode<BniNetwork> I2 = new BniNetNode<>();
        NetNode<BniNetwork> I3 = new BniNetNode<>();
        NetNode<BniNetwork> I4 = new BniNetNode<>();


        root.adoptChild(I1, 2.0);
        root.adoptChild(I2, 3.0);
        I1.adoptChild(I2, 2.0);

        I1.adoptChild(I3, 4.0);
        I2.adoptChild(I4, 4.0);
        I3.adoptChild(I4, 2.0);

        I3.adoptChild(b, 4.0);
        I4.adoptChild(a, 3.0);

        return testNet;
    }

    public static void main(String[] args) {
        Network<BniNetwork> net1 = getNetwork1();

        System.out.println("Topology Trees: ");
        for (NetworkTree<BniNetwork> tree : Networks.getTrees(net1)) {
            System.out.println(tree.makeTree());
        }

        System.out.println();

        System.out.println("Length Trees: ");
        for (STITree tree : Networks.getExtendedTrees(net1)) {
            System.out.println(tree);
        }

        System.out.println();

        System.out.println("Net1: " + net1);

        System.out.println(Arrays.deepToString(AveragePathDistance.computeMatrix(net1)));
    }
}
