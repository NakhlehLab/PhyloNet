package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal.tests;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.junit.Test;

import java.util.Random;

/**
 * Created by wendingqiao on 10/20/14.
 */
public class NetworkBranchScalerTest {
    @Test
    public void test() throws Exception
    {
        Random random = new Random(294345L);

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

        k.adoptChild(d, 2.0);
        j.adoptChild(a, 2.0);
        k.adoptChild(x, 0.0);
        j.adoptChild(x, 0.0);
        x.adoptChild(l, 1.0);
        l.adoptChild(b, 1.0);
        l.adoptChild(c, 1.0);

        x.setParentProbability(j, .3);
        x.setParentProbability(k,.7);


        // test undo()
        for(int i = 0; i < 1000; i++) {
            NetInheritScaler op = new NetInheritScaler<>(network, random);
            System.out.println(op.operate());
            System.out.println(network.toString());
            op.undo();
            System.out.println(network.toString());
        }
        // test operation()
        for(int i = 0; i < 1000; i++) {
            NetInheritScaler op = new NetInheritScaler<>(network, random);
            System.out.println(op.operate());
            System.out.println(network.toString());
        }


        // test undo()

        System.out.println(network.toString());
        for(int i = 0; i < 1000; i++) {
            NetworkOperation<Double> op = new NetBranchScaler<Double>(network, random);
            System.out.println(op.operate());
            //System.out.println(network.toString());
            op.undo();
            //System.out.println(network.toString());
        }
        System.out.println(network.toString());


        // test operation()
        for(int i = 0; i < 1000; i++) {
            NetBranchScaler op = new NetBranchScaler<>(network, random);
            System.out.println(op.operate());
            System.out.println(network.toString());
        }

    }
}

