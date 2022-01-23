package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.junit.Assert;
import org.junit.Test;


public class TreeBasedDissimilarityTest {

    /**
     * Test if the score is 0 if the networks are the same.
     */
    @Test
    public void testNormalizedTreeDistanceZeroDifference() {
        String richNewick = "((((((C:1.0)#H1:1.0::0.5,B:2.0):2.0)#H2:1.0::0.5,A:5.0):3.0,((D:3.0,#H1:2.0::0.5):3.0,#H2:2.0::0.5):2.0):2.0,E:10.0);";

        Network<BniNetwork> network1 = Networks.readNetwork(richNewick);
        Network<BniNetwork> network2 = Networks.readNetwork(richNewick);

        TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity<>(network1, network2);

        Assert.assertEquals(0.0, treeBasedDissimilarity.computeNormalizedTreeDistance(), 1e-5);
    }

    /**
     * Test if the score is 0 if the networks are the same.
     */
    @Test
    public void testRootedBranchScoreZeroDifference() {
        String richNewick = "((((((C:1.0)#H1:1.0::0.5,B:2.0):2.0)#H2:1.0::0.5,A:5.0):3.0,((D:3.0,#H1:2.0::0.5):3.0,#H2:2.0::0.5):2.0):2.0,E:10.0);";

        Network<BniNetwork> network1 = Networks.readNetwork(richNewick);
        Network<BniNetwork> network2 = Networks.readNetwork(richNewick);

        TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity<>(network1, network2);

        Assert.assertEquals(0.0, treeBasedDissimilarity.computeRootedBranchScore(), 1e-5);
    }

    /**
     * Checks whether our network measure behaves exactly same as the rooted branch score if there are no reticulations.
     */
    @Test
    public void testRootedBranchScoreConsistency() {
        String newick1 = "(E:10,((B:3,(F:1,G:1):2):4,(D:6,(H:4,(C:2,A:2):2):2):1):3);";
        String newick2 = "(E:20,((B:6,(F:2,G:2):4):8,(D:12,(H:8,(C:4,A:4):4):4):2):6);";

        Network<BniNetwork> network1 = Networks.readNetwork(newick1);
        Network<BniNetwork> network2 = Networks.readNetwork(newick2);

        Tree tree1 = Trees.readTree(newick1);
        Tree tree2 = Trees.readTree(newick2);

        TreeBasedDissimilarity<BniNetwork> treeBasedDissimilarity = new TreeBasedDissimilarity<>(network1, network2);
        double dissimilarity1 = treeBasedDissimilarity.computeRootedBranchScore();
        double dissimilarity2 = TreeBasedDissimilarity.computeRootedBranchScore(tree1, 1, tree2, 1);

        Assert.assertEquals(dissimilarity1, dissimilarity2, 1e-5);
    }
}