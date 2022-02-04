package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.junit.Assert;
import org.junit.Test;


public class RootedNetworkBranchScoreTest {

    /**
     * Test if the score is 0 if the networks are the same.
     */
    @Test
    public void testZeroDifference() {
        String richNewick = "((((((C:1.0)#H1:1.0::0.5,B:2.0):2.0)#H2:1.0::0.5,A:5.0):3.0,((D:3.0,#H1:2.0::0.5):3.0,#H2:2.0::0.5):2.0):2.0,E:10.0);";

        Network<BniNetwork> network1 = Networks.readNetwork(richNewick);
        Network<BniNetwork> network2 = Networks.readNetwork(richNewick);

        RootedNetworkBranchScore<BniNetwork> RNBS = new RootedNetworkBranchScore<>(network1, network2);

        Assert.assertEquals(0.0, RNBS.compute(), 1e-5);
        Assert.assertEquals(0.0, RNBS.computeNormalized(), 1e-5);
    }
}