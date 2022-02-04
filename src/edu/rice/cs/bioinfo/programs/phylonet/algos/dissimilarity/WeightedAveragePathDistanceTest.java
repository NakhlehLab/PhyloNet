package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.junit.Assert;
import org.junit.Test;

public class WeightedAveragePathDistanceTest {
    @Test
    public void testZeroDifference() {
        String richNewick1 = "(E:10.0,((((B:2.0,(C:1.0)I7#H2:1.0::0.5)I6:2.0)I4#H1:2.0::0.5,(I7#H2:2.0::0.5,D:3.0)I5:3.0)I3:2.0,(A:5.0,I4#H1:1.0::0.5)I2:3.0)I1:2.0)I0;";
        String richNewick2 = "(E:10.0,((((B:2.0,(C:1.0)I7#H2:1.0::0.5)I6:2.0)I4#H1:2.0::0.5,(I7#H2:2.0::0.5,D:3.0)I5:3.0)I3:2.0,(A:5.0,I4#H1:1.0::0.5)I2:3.0)I1:2.0)I0;";

        Network<BniNetwork> network1 = Networks.readNetwork(richNewick1);
        Network<BniNetwork> network2 = Networks.readNetwork(richNewick2);

        WeightedAveragePathDistance<BniNetwork> WAPD = new WeightedAveragePathDistance<>(network1, network2);

        Assert.assertEquals(0.0, WAPD.compute(), 1e-5);
        Assert.assertEquals(0.0, WAPD.computeNormalized(), 1e-5);
    }
}
