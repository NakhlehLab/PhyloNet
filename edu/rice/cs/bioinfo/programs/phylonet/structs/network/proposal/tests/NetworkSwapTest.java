package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal.tests;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal.NetworkSwap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.junit.Test;

import java.util.Random;

/**
 * Created by wendingqiao on 10/18/14.
 */
public class NetworkSwapTest {
    @Test
    public void test() throws Exception
    {
        Random random = new Random(475905L);

        BniNetwork<Integer> network = new BniNetwork<Integer>();
        network.createRoot("R");

        NetNode<Integer> r = network.getRoot();
        NetNode<Integer> j = new BniNetNode<Integer>("J", null);
        NetNode<Integer> k = new BniNetNode<Integer>("K", null);
        NetNode<Integer> a = new BniNetNode<Integer>("A", null);
        NetNode<Integer> x = new BniNetNode<Integer>("X", null);
        NetNode<Integer> d = new BniNetNode<Integer>("D", null);
        NetNode<Integer> l = new BniNetNode<Integer>("L", null);
        NetNode<Integer> b = new BniNetNode<Integer>("B", null);
        NetNode<Integer> c = new BniNetNode<Integer>("C", null);

        r.adoptChild(j, 1);
        r.adoptChild(k, 1);

        k.adoptChild(d, 2);
        j.adoptChild(a, 2);
        k.adoptChild(x, 0);
        j.adoptChild(x, 0);
        x.adoptChild(l, 1);
        l.adoptChild(b, 1);
        l.adoptChild(c, 1);

        x.setParentProbability(j, .3);
        x.setParentProbability(k,.7);

        // test undo()
        System.out.println(network.toString());
        for(int i = 0; i < 1000; i++) {
            NetworkSwap op = new NetworkSwap(network, random);
            op.operate();
            op.undo();
        }
        System.out.println(network.toString());
        // use dendroscope to test the equality of the two figures.

        // test operation()
        for(int i = 0; i < 1000; i++) {
            NetworkSwap op = new NetworkSwap(network, random);
            System.out.println(op.operate());
            assert(!Networks.hasCycle(network));
            System.out.println(network.toString());
        }


        /*
        // test cycle
        BniNetwork<Integer> cycle = new BniNetwork<Integer>();
        cycle.createRoot("A");
        NetNode<Integer> aa = cycle.getRoot();
        NetNode<Integer> bb = new BniNetNode<Integer>("B", null);
        NetNode<Integer> cc = new BniNetNode<Integer>("C", null);
        NetNode<Integer> dd = new BniNetNode<Integer>("D", null);
        NetNode<Integer> ee = new BniNetNode<Integer>("E", null);
        aa.adoptChild(bb, 1.0);
        aa.adoptChild(dd, 1.0);
        bb.adoptChild(cc, 1.0);
        bb.adoptChild(ee, 1.0);
        dd.adoptChild(cc, 1.0);
        ee.adoptChild(cc, 1.0);
        cc.adoptChild(ee, 1.0);
        assert(Networks.hasCycle(cycle));
        */
    }
}
