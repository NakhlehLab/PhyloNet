package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal.tests;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal.NetRetAdd;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal.NetRetDel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal.NetworkOperation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal.NetworkSwap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.junit.Test;

import java.util.Random;

/**
 * Created by wendingqiao on 10/21/14.
 */
public class NetworkReticulationTest {
    @Test
    public void test() throws Exception {

        Random random = new Random(559872L);

        BniNetwork<Double> network = new BniNetwork<Double>();
        network.createRoot("R");
        NetNode<Double> r = network.getRoot();
        NetNode<Double> j = new BniNetNode<Double>("J", null);
        NetNode<Double> k = new BniNetNode<Double>("K", null);
        NetNode<Double> a = new BniNetNode<Double>("A", null);
        NetNode<Double> x = new BniNetNode<Double>("X", null);
        NetNode<Double> d = new BniNetNode<Double>("D", null);
        NetNode<Double> l = new BniNetNode<Double>("L", null);
        NetNode<Double> b = new BniNetNode<Double>("B", null);
        NetNode<Double> c = new BniNetNode<Double>("C", null);

        r.adoptChild(j, 1.0);
        r.adoptChild(k, 1.0);
        k.adoptChild(d, 1.0);
        j.adoptChild(a, 1.0);
        k.adoptChild(x, 1.0);
        j.adoptChild(x, 1.0);
        x.adoptChild(l, 1.0);
        l.adoptChild(b, 1.0);
        l.adoptChild(c, 1.0);

        x.setParentProbability(j, 0.5);
        x.setParentProbability(k, 0.5);





        BniNetwork<Double> net1 = new BniNetwork<Double>();
        net1.createRoot("R");
        NetNode<Double> R = net1.getRoot();
        NetNode<Double> J = new BniNetNode<Double>("J", null);
        NetNode<Double> K = new BniNetNode<Double>("K", null);
        NetNode<Double> A = new BniNetNode<Double>("A", null);
        NetNode<Double> X = new BniNetNode<Double>("X", null);
        NetNode<Double> D = new BniNetNode<Double>("D", null);
        NetNode<Double> L = new BniNetNode<Double>("L", null);
        NetNode<Double> B = new BniNetNode<Double>("B", null);
        NetNode<Double> C = new BniNetNode<Double>("C", null);

        R.adoptChild(J, 1.0);
        R.adoptChild(K, 1.0);
        K.adoptChild(D, 1.0);
        J.adoptChild(A, 1.0);
        K.adoptChild(X, 1.0);
        J.adoptChild(X, 1.0);
        X.adoptChild(L, 1.0);
        L.adoptChild(B, 1.0);
        L.adoptChild(C, 1.0);

        X.setParentProbability(J, 0.5);
        X.setParentProbability(K, 0.5);


        // Network reticulation node addition
        // test undo()
        System.out.println(network.toString());
        for (int i = 0; i < 0; i++) {
            NetworkOperation op = new NetRetAdd<>(network, random);
            System.out.println(op.operate());
            //System.out.println(network.toString());
            op.undo();
            //System.out.println(network.toString());
        }
        System.out.println(network.toString());
        // test operation()
        for (int i = 0; i < 100; i++) {
            //System.out.println(network.toString());
            NetworkOperation op = new NetRetAdd<>(network, random);
            System.out.println(op.operate());
            assert (!Networks.hasCycle(network));
            //op = new NetRetDel<>(network);
            //System.out.println(op.operate());
            System.out.println(Networks.hasTheSameTopology(network, net1));
            System.out.println(network.toString());
        }

        // Network reticulation node deletion

        // test undo()
        System.out.println(network.toString());
        for (int i = 0; i < 100; i++) {
            //System.out.println(network.toString());
            NetworkOperation op = new NetRetDel<>(network, random);
            System.out.println(op.operate());
            //System.out.println(network.toString());
            op.undo();
            //System.out.println(network.toString());
        }
        System.out.println(network.toString());
        // use dendroscope to test the equality of the two figures.

        // test operate()
        for (int i = 0; i < 100; i++) {
            System.out.println(network.toString());
            NetworkOperation op = new NetRetDel<>(network, random);
            System.out.println(op.operate());
            System.out.println(Networks.hasTheSameTopology(network, net1));
            //System.out.println(network.toString());
        }
    }
}
