package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
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

public class TreeBasedDissimilarity<T> {
    private Network<T> network1;
    private Network<T> network2;

    private Iterable<NetworkTree<T>> networkTrees1;
    private Iterable<NetworkTree<T>> networkTrees2;

    private List<Tree> trees1;
    private List<Tree> trees2;

    /**
     * Constructor for Tree-Based Dissimilarity metric.
     * @param network1 First network.
     * @param network2 Second network.
     */
    public TreeBasedDissimilarity(Network<T> network1, Network<T> network2){
        this.network1 = network1;
        this.network2 = network2;
        this.networkTrees1 = Networks.getTrees(this.network1);
        this.networkTrees2 = Networks.getTrees(this.network2);

        this.trees1 = new ArrayList<>();
        this.trees2 = new ArrayList<>();
        this.networkTrees1.iterator().forEachRemaining(i-> this.trees1.add(i.makeTree()));
        this.networkTrees2.iterator().forEachRemaining(i-> this.trees2.add(i.makeTree()));
    }

    /**
     * Computes the distance in between two networks using tree-based dissimilarity metric.
     *
     * @return A dissimilarity measure.
     */
    public double computeNormalizedTreeDistance() {
        // Initialize a bipartite graph with nodes corresponding to trees in the two networks.
        BipartiteGraph BG = new BipartiteGraph(trees1.size(), trees2.size());

        for (int l = 0; l < trees1.size(); l++) {
            for (int r = 0; r < trees2.size(); r++) {
                if (!Trees.haveSameRootedTopology(trees1.get(l), trees2.get(r))) {
                    continue;
                }
                double weight = computeNormalizedTreeDistance(trees1.get(l), 1.0, trees2.get(r), 1.0);
                BG.addEdge(l, r, weight);
            }
        }

        return BG.getMinEdgeCoverWeight();
    }

    /**
     * Computes the distance in between two networks using tree-based dissimilarity metric.
     *
     * @return A dissimilarity measure.
     */
    public double computeRootedBranchScore() {
        // Initialize a bipartite graph with nodes corresponding to trees in the two networks.
        BipartiteGraph BG = new BipartiteGraph(trees1.size(), trees2.size());

        for (int l = 0; l < trees1.size(); l++) {
            for (int r = 0; r < trees2.size(); r++) {
                double weight = computeRootedBranchScore(trees1.get(l), 1.0, trees2.get(r), 1.0);
                BG.addEdge(l, r, weight);
            }
        }

        return BG.getMinEdgeCoverWeight();
    }

    /**
     * Compute normalized tree distance for trees tree1 and tree2. The method is described in:
     * Zheng, Yichen et al. A scale-free method for testing the proportionality of branch lengths between two
     * phylogenetic trees. arXiv: Quantitative Methods (2015): n. pag.
     *
     * @param tree1 First tree.
     * @param scale1 Scale factor for branch lengths of tree1.
     * @param tree2 Second tree.
     * @param scale2 Scale factor for branch lengths of tree2.
     * @return Normalized tree distance.
     */
     static double computeNormalizedTreeDistance(Tree tree1, double scale1, Tree tree2, double scale2) {
        if (!Trees.leafSetsAgree(tree1, tree2)) {
            throw new RuntimeException("Trees must have identical leaf sets");
        }

        if (!Trees.haveSameRootedTopology(tree1, tree2)) {
            throw new RuntimeException("Trees must have identical topology");
        }

        List<Double> tree1BranchLengths = new ArrayList<>();
        double tree1BranchLengthSum = 0;
        for (TNode node : tree1.postTraverse()) {
            if (!node.isRoot()) {
                tree1BranchLengthSum += node.getParentDistance() * scale1;
                tree1BranchLengths.add(node.getParentDistance() * scale1);
            }
        }

        List<Double> tree2BranchLengths = new ArrayList<>();
        double tree2BranchLengthSum = 0;
        for (TNode node : tree2.postTraverse()) {
            if (!node.isRoot()) {
                tree2BranchLengthSum += node.getParentDistance() * scale2;
                tree2BranchLengths.add(node.getParentDistance() * scale2);
            }
        }

        if (tree1BranchLengths.size() != tree2BranchLengths.size()) {
            throw new RuntimeException("Trees must have same number of branches");
        }

        if (Double.isInfinite(tree1BranchLengthSum) || Double.isNaN(tree1BranchLengthSum))
            throw new RuntimeException("Branch length cannot be infinite or undefined");

        if (Double.isInfinite(tree2BranchLengthSum) || Double.isNaN(tree2BranchLengthSum))
            throw new RuntimeException("Branch length cannot be infinite or undefined");

        double sum = 0;
        for (int i = 0; i < tree1BranchLengths.size(); i++) {
            sum += Math.abs(tree1BranchLengths.get(i) / tree1BranchLengthSum - tree2BranchLengths.get(i) / tree2BranchLengthSum);
        }

        return sum / 2.0;
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
     static double computeRootedBranchScore(Tree tree1, double scale1, Tree tree2, double scale2) {
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