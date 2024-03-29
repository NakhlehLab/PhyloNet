package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.BipartiteGraph;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class RootedNetworkBranchScore {

    /**
     * Compute rNBS.
     *
     * @return rNBS of two networks.
     */
    public static <T> double compute(Network<T> network1, Network<T> network2) {
        List<Tree> trees1 = new ArrayList<>();
        List<Tree> trees2 = new ArrayList<>();
        Networks.getExtendedTrees(network1).iterator().forEachRemaining(trees1::add);
        Networks.getExtendedTrees(network2).iterator().forEachRemaining(trees2::add);

        // Initialize a bipartite graph with nodes corresponding to trees in the two networks.
        BipartiteGraph BG = new BipartiteGraph(trees1.size(), trees2.size());

        for (int l = 0; l < trees1.size(); l++) {
            for (int r = 0; r < trees2.size(); r++) {
                double weight = computeRootedBranchScore(trees1.get(l), 1.0, trees2.get(r), 1.0);
                BG.addEdge(l, r, weight);
            }
        }

        return BG.getMinEdgeCoverWeight() / BG.getMinEdgeCoverSize();
    }

    /**
     * Compute rNBS.
     *
     * @return rNBS of two networks.
     */
    public static double compute(List<Tree> trees1, List<Tree> trees2) {
        // Initialize a bipartite graph with nodes corresponding to trees in the two networks.
        BipartiteGraph BG = new BipartiteGraph(trees1.size(), trees2.size());

        for (int l = 0; l < trees1.size(); l++) {
            for (int r = 0; r < trees2.size(); r++) {
                double weight = computeRootedBranchScore(trees1.get(l), 1.0, trees2.get(r), 1.0);
                BG.addEdge(l, r, weight);
            }
        }

        return BG.getMinEdgeCoverWeight() / BG.getMinEdgeCoverSize();
    }

    /**
     * Compute normalized rNBS.
     *
     * @return Normalized rNBs of two networks.
     */
    public static <T> double computeNormalized(Network<T> network1, Network<T> network2) {
        List<Tree> trees1 = new ArrayList<>();
        List<Tree> trees2 = new ArrayList<>();
        Networks.getExtendedTrees(network1).iterator().forEachRemaining(trees1::add);
        Networks.getExtendedTrees(network2).iterator().forEachRemaining(trees2::add);

        // Initialize a bipartite graph with nodes corresponding to trees in the two networks.
        BipartiteGraph BG = new BipartiteGraph(trees1.size(), trees2.size());

        for (int l = 0; l < trees1.size(); l++) {
            for (int r = 0; r < trees2.size(); r++) {
                double weight = computeRootedBranchScore(trees1.get(l), 1.0, trees2.get(r), 1.0);
                weight /= Math.min(Trees.getTotalBranchLength(trees1.get(l)), Trees.getTotalBranchLength(trees2.get(r)));
                BG.addEdge(l, r, weight);
            }
        }

        return BG.getMinEdgeCoverWeight() / BG.getMinEdgeCoverSize();
    }

    /**
     * Compute normalized rNBS.
     *
     * @return Normalized rNBs of two networks.
     */
    public static double computeNormalized(List<Tree> trees1, List<Tree> trees2) {
        // Initialize a bipartite graph with nodes corresponding to trees in the two networks.
        BipartiteGraph BG = new BipartiteGraph(trees1.size(), trees2.size());

        for (int l = 0; l < trees1.size(); l++) {
            for (int r = 0; r < trees2.size(); r++) {
                double weight = computeRootedBranchScore(trees1.get(l), 1.0, trees2.get(r), 1.0);
                weight /= Math.min(Trees.getTotalBranchLength(trees1.get(l)), Trees.getTotalBranchLength(trees2.get(r)));
                BG.addEdge(l, r, weight);
            }
        }

        return BG.getMinEdgeCoverWeight() / BG.getMinEdgeCoverSize();
    }

    /**
     * Compute rooted branch score for trees tree1 and tree2. The method is proposed in:
     * M K Kuhner, J Felsenstein, A simulation comparison of phylogeny algorithms under equal and unequal evolutionary
     * rates., Molecular Biology and Evolution, Volume 11, Issue 3, May 1994, Pages 459–468,
     * https://doi.org/10.1093/oxfordjournals.molbev.a040126
     *
     * and described in:
     * Joseph Heled, Alexei J. Drummond, Bayesian Inference of Species Trees from Multilocus Data, Molecular Biology and
     * Evolution, Volume 27, Issue 3, March 2010, Pages 570–580, https://doi.org/10.1093/molbev/msp274
     *
     * @param tree1 First tree.
     * @param scale1 Scale factor for branch lengths of tree1.
     * @param tree2 Second tree.
     * @param scale2 Scale factor for branch lengths of tree2.
     * @return Rooted branch score.
     */
     private static double computeRootedBranchScore(Tree tree1, double scale1, Tree tree2, double scale2) {
        if (!Trees.leafSetsAgree(tree1, tree2)) {
            throw new RuntimeException("Trees must have identical leaf sets");
        }
        String[] taxa = tree1.getLeaves();
        Map<STITreeCluster, TNode> clusters1 = ((STITree)tree1).getClusters(taxa);
        Map<STITreeCluster, TNode> clusters2 = ((STITree)tree2).getClusters(taxa);

        double sum = 0;

        for (STITreeCluster cl1: clusters1.keySet()) {
            double dist;

            if (clusters2.containsKey(cl1)) {
                dist = clusters1.get(cl1).getParentDistance() * scale1 - clusters2.get(cl1).getParentDistance() * scale2;
            } else {
                dist = clusters1.get(cl1).getParentDistance() * scale1;
            }

            if (Double.isInfinite(dist) || Double.isNaN(dist))
                throw new RuntimeException("Branch length cannot be infinite or undefined");

            sum += dist * dist;
        }

        for (STITreeCluster cl2: clusters2.keySet()) {
            double dist;

            if (clusters1.containsKey(cl2)) {
                continue; // Has been processed.
            }

            dist = clusters2.get(cl2).getParentDistance() * scale2;

            if (Double.isInfinite(dist) || Double.isNaN(dist))
                throw new RuntimeException("Branch length cannot be infinite or undefined");

            sum += dist * dist;
        }

        return Math.sqrt(sum);
    }
}