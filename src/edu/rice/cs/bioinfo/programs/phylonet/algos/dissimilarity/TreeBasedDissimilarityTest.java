package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.Random;

public class TreeBasedDissimilarityTest {
    private Random random;

    @Before
    public void init() {
        this.random = new Random();
    }

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

    @Test
    public void testRootedBranchScoreConsistency() {
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

        String[] letters = new String[8];
        for (int i = 0; i < 8; i++) {
            letters[i] = String.valueOf(alphabet.charAt(i));
        }

        List<Tree> newTrees = Trees.generateAllBinaryTrees(letters);
        Tree newTree = newTrees.get(this.random.nextInt(newTrees.size()));

        System.out.println(newTree);

        String richNewick1 = "((((((C:1.0)#H1:1.0::0.5,B:2.0):2.0)#H2:1.0::0.5,A:5.0):3.0,((D:3.0,#H1:2.0::0.5):3.0,#H2:2.0::0.5):2.0):2.0,E:10.0);";

        Network<BniNetwork> network1 = Networks.readNetwork(richNewick1);
    }
}