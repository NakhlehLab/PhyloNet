package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.junit.Assert;
import org.junit.Test;

public class TreeBasedDissimilarityTest {
    @Test
    public void testNormalizedTreeDistanceZeroDifference() {
        String richNewick1 = "((((((C:1.0)#H1:1.0::0.5,B:2.0):2.0)#H2:1.0::0.5,A:5.0):3.0,((D:3.0,#H1:2.0::0.5):3.0,#H2:2.0::0.5):2.0):2.0,E:10.0);";
        String richNewick2 = "((((((C:1.0)#H1:1.0::0.5,B:2.0):2.0)#H2:1.0::0.5,A:5.0):3.0,((D:3.0,#H1:2.0::0.5):3.0,#H2:2.0::0.5):2.0):2.0,E:10.0);";

        Network<BniNetwork> network1 = Networks.readNetwork(richNewick1);
        Network<BniNetwork> network2 = Networks.readNetwork(richNewick2);

        TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity<>(network1, network2);

        Assert.assertEquals(0.0, treeBasedDissimilarity.computeNormalizedTreeDistance(), 1e-5);
    }

    @Test
    public void testRootedBranchScoreZeroDifference() {
        String richNewick1 = "((((((C:1.0)#H1:1.0::0.5,B:2.0):2.0)#H2:1.0::0.5,A:5.0):3.0,((D:3.0,#H1:2.0::0.5):3.0,#H2:2.0::0.5):2.0):2.0,E:10.0);";
        String richNewick2 = "((((((C:1.0)#H1:1.0::0.5,B:2.0):2.0)#H2:1.0::0.5,A:5.0):3.0,((D:3.0,#H1:2.0::0.5):3.0,#H2:2.0::0.5):2.0):2.0,E:10.0);";

        Network<BniNetwork> network1 = Networks.readNetwork(richNewick1);
        Network<BniNetwork> network2 = Networks.readNetwork(richNewick2);

        TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity<>(network1, network2);

        Assert.assertEquals(0.0, treeBasedDissimilarity.computeRootedBranchScore(), 1e-5);
    }
}