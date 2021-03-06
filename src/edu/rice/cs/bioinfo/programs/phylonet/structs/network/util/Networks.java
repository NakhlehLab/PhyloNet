/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.util;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.fitchpars.ParsimonyCalculator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTreeEnumerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTripartition;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.sequence.model.SequenceAlignment;
import edu.rice.cs.bioinfo.programs.phylonet.structs.sequence.model.SequenceException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.ByteArrayInputStream;
import java.math.BigDecimal;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/31/11
 * Time: 3:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class Networks
{
    /**
     * This function names internal nodes so that they all have non-empty names.
     * It assumes that there are currently no duplicate names in the network.
     */
    public static <T> void autoLabelNodes(Network<T> net)
    {
        Map<String, NetNode<T>> map = new Hashtable<String, NetNode<T>>();
        for (NetNode<T> node : net.bfs()) {
            if (!node.getName().equals(NetNode.NO_NAME)) {
                map.put(node.getName(), node);
            }
        }

        int number = 0;
        for (NetNode<T> node : net.bfs()) {
            if (node.getName().equals(NetNode.NO_NAME)) {
                String name;
                do {
                    name = NAME_PREFIX + (number++);
                } while (map.get(name) != null);
                node.setName(name);
            }
        }
    }


    /**
     * This function removes all binary nodes from the network.
     */
    public static <T> void removeBinaryNodes(Network<T> net) {
        boolean update;
        do{
            update = false;
            for (NetNode<T> node : net.bfs()) {
                if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                    NetNode<T> parent = node.getParents().iterator().next();    // Node's only parent.
                    NetNode<T> child = node.getChildren().iterator().next();    // Node's only child.
                    double distance = node.getParentDistance(parent) + child.getParentDistance(node);
                    double support = Math.max(node.getParentSupport(parent), child.getParentSupport(node));
                    double inheritanceProb = child.getParentProbability(node);
                    parent.removeChild(node);
                    node.removeChild(child);
                    parent.adoptChild(child, distance);
                    child.setParentProbability(parent, inheritanceProb);
                    child.setParentSupport(parent, support);
                    update = true;
                }
                if(node.isRoot() && node.getOutdeg()==1){
                    NetNode newRoot = node.getChildren().iterator().next();
                    double root_popsize = node.getRootPopSize();
                    node.removeChild(newRoot);
                    net.resetRoot(newRoot);
                    newRoot.setRootPopSize(root_popsize);
                }
            }
        }while(update);
    }


    /**
     * This function creates a dummy network, which contains just the root and a set of leaves
     * directly connected to the root. This function is just for creating a NetworkCluster from
     * a set of leaves when the actual network is unknown.
     */
    public static <T> Network<T> createDummyNetwork(List<String> leaf_names)
    {
        BniNetNode<T> root = new BniNetNode<T>();
        for (String name : leaf_names) {
            BniNetNode<T> leaf = new BniNetNode<T>();
            leaf.setName(name);
            root.adoptChild(leaf, NetNode.NO_DISTANCE);
        }

        return new BniNetwork<T>(root);
    }


    /**
     * This function returns a list of trees decomposed from the network.
     *
     * @param net: The network from which trees are generated.
     */
    public static <T> Iterable<NetworkTree<T>> getTrees(Network<T> net)
    {
        List<NetworkTree<T>> trees = new LinkedList<NetworkTree<T>>();

        if (!net.isEmpty()) {
            NetworkTreeEnumerator<T> enumerator = new NetworkTreeEnumerator<T>(net);
            for (NetworkTree<T> nt : enumerator) {
                if (!trees.contains(nt)) {
                    trees.add(nt);
                }
            }
        }

        return trees;
    }

    /**
     * This function returns a list of trees decomposed from the network that display unique lengths.
     * @param net: The network from which trees are generated.
     * @return Iterable object over trees.
     */
    public static <T> Iterable<Tree> getExtendedTrees(Network<T> net)
    {
        List<Tree> trees = new LinkedList<>();

        if (!net.isEmpty()) {
            NetworkTreeEnumerator<T> enumerator = new NetworkTreeEnumerator<T>(net);
            for (NetworkTree<T> nt : enumerator) {
                Tree t = nt.makeTree();
                t.getRoot().setParentDistance(NetNode.NO_DISTANCE); // TODO ALP: Investigate why..
                if (!trees.contains(t)) {
                    trees.add(t);
                }
            }
        }

        return trees;
    }

    /**
     * This function returns an iterable list of network clusters.
     *
     * @param net: The network we want to find all clusters.
     */
    public static <T> Iterable<NetworkCluster<T>> getSoftwiredClusters(Network<T> net)
    {
        List<NetworkCluster<T>> clusters = new LinkedList<NetworkCluster<T>>();

        for (NetworkTree<T> nt : getTrees(net)) {
            for (NetworkCluster<T> nc : nt.getClusters()) {
                if (!clusters.contains(nc)) {
                    clusters.add(nc);
                }
            }
        }

        return clusters;
    }

    /**
     * This function returns an iterable list of network clusters.
     *
     * @param net: The network we want to find all clusters.
     */
    public static <T> Iterable<Set<String>> getHardwiredClusters(Network<T> net)
    {

        Map<NetNode, Set<String>> node2cluster = new HashMap<NetNode, Set<String>>();
        for (NetNode<T> node: Networks.postTraversal(net)) {
            if(node.isRoot()){
                break;
            }
            Set<String> cluster = new HashSet<String>();
            if(node.isLeaf()){
                cluster.add(node.getName());
            }
            else{
                for(NetNode<T> child: node.getChildren()){
                    cluster.addAll(node2cluster.get(child));
                }
            }
            node2cluster.put(node, cluster);
        }
        Set<Set<String>> clusters = new HashSet<Set<String>>();
        for(Set<String> cluster: node2cluster.values()){
            if(cluster.size()>1){
                clusters.add(cluster);
            }
        }

        return clusters;
    }


    /**
     * This function returns an iterable list of network tripartitions.
     */
    public static <T> Iterable<NetworkTripartition<T>> getTripartitions(Network<T> net)
    {
        List<NetworkTripartition<T>> tripartitions = new LinkedList<NetworkTripartition<T>>();

        for (NetNode<T> node : net.bfs()) {
            if (!node.isRoot() && !node.isLeaf()) {
                NetworkTripartition<T> ntp = new NetworkTripartition<T>(net, node);

                if (!tripartitions.contains(ntp)) {
                    tripartitions.add(ntp);
                }
            }
        }

        return tripartitions;
    }

    /**
     * This function computes the softwired-cluster distance between two networks.
     *
     * @param net1, net2: The two networks to be compared.
     *
     * @return [fasle-negative, false-positive, average]
     */
    public static <T> double[] computeSoftwiredClusterDistance(Network<T> net1, Network<T> net2)
    {
        List<NetworkCluster<T>> clusters1 = new LinkedList<NetworkCluster<T>>();
        List<NetworkCluster<T>> clusters2 = new LinkedList<NetworkCluster<T>>();

        for (NetworkCluster<T> nc : getSoftwiredClusters(net1)) {
            clusters1.add(nc);
        }
        for (NetworkCluster<T> nc : getSoftwiredClusters(net2)) {
            clusters2.add(nc);
        }

        return computeSoftwiredClusterDistance(clusters1, clusters2);
    }


    /**
     * This function computes the softwired-cluster distance given two lists of softwired clusters
     *
     * @param clusters1, clusters2: The two lists of softwired clusters to be compared.
     *
     * @return [fasle-negative, false-positive, average]
     */
    public static <T> double[] computeSoftwiredClusterDistance(List<NetworkCluster<T>> clusters1, List<NetworkCluster<T>> clusters2)
    {
        double fn = (clusters1.size() == 0) ? 0.0 : (double) computeSoftwiredClusterDiff(clusters1, clusters2) / clusters1.size();
        double fp = (clusters2.size() == 0) ? 0.0 : (double) computeSoftwiredClusterDiff(clusters2, clusters1) / clusters2.size();
        double avg = (fn + fp) / 2;

        return new double[] {fn, fp, avg};
    }




    /**
     * This function computes the number of softwired clusters in <code>models</code> that are not in <code>refs</code>.
     */
    private static <T> int computeSoftwiredClusterDiff(List<NetworkCluster<T>> models, List<NetworkCluster<T>> refs)
    {
        int diff = 0;
        for (NetworkCluster<T> nc : models) {
            if (!refs.contains(nc)) {
                diff++;
            }
        }

        return diff;
    }


    /**
     * This function compute the hardwired-cluster distance between two networks.
     *
     * @param net1, net2: The two networks to be compared.
     *
     * @return [fasle-negative, false-positive, average]
     */
    public static <T> double[] computeHardwiredClusterDistance(Network<T> net1, Network<T> net2)
    {
        List<Set<String>> clusters1 = new LinkedList<Set<String>>();
        List<Set<String>> clusters2 = new LinkedList<Set<String>>();

        for (Set<String> cl : getHardwiredClusters(net1)) {
            clusters1.add(cl);
        }
        for (Set<String> cl : getHardwiredClusters(net2)) {
            clusters2.add(cl);
        }

        return computeHardwiredClusterDistance(clusters1, clusters2);
    }



    /**
     * This function computes the hardwired-cluster distance given two lists of softwired clusters
     *
     * @param clusters1, clusters2: The two lists of softwired clusters to be compared.
     *
     * @return [fasle-negative, false-positive, average]
     */
    public static <T> double[] computeHardwiredClusterDistance(List<Set<String>> clusters1, List<Set<String>> clusters2)
    {
        double fn = (clusters1.size() == 0) ? 0.0 : (double) computeHardwiredClusterDiff(clusters1, clusters2) / clusters1.size();
        double fp = (clusters2.size() == 0) ? 0.0 : (double) computeHardwiredClusterDiff(clusters2, clusters1) / clusters2.size();
        double avg = (fn + fp) / 2;

        return new double[] {fn, fp, avg};
    }



    /**
     * This function computes the number of hardwired clusters in <code>models</code> that are not in <code>refs</code>.
     */
    private static <T> int computeHardwiredClusterDiff(List<Set<String>> models, List<Set<String>> refs)
    {
        int diff = 0;
        for (Set<String> nc : models) {
            if (!refs.contains(nc)) {
                diff++;
            }
        }

        return diff;
    }


    /**
     * This function computes the tripartition-based distance between two networks.
     *
     * @param net1, net2: Two networks to be compared.
     *
     * @return [fasle-negative, false-positive, average]
     */
    public static <T> double[] computeTripartitionDistance(Network<T> net1, Network<T> net2)
    {
        List<NetworkTripartition<T>> partitions1 = new LinkedList<NetworkTripartition<T>>();
        List<NetworkTripartition<T>> partitions2 = new LinkedList<NetworkTripartition<T>>();

        for (NetworkTripartition<T> ntp : getTripartitions(net1)) {
            partitions1.add(ntp);
        }
        for (NetworkTripartition<T> ntp : getTripartitions(net2)) {
            partitions2.add(ntp);
        }

        return computeTripartitionDistance(partitions1, partitions2);
    }



    /**
     * This function computes the tripartition-based distance given two lists of partitions
     *
     * @param partitions1, partitions2: Two lists of partitions to be compared.
     *
     * @return [fasle-negative, false-positive, average]
     */
    public static <T> double[] computeTripartitionDistance(List<NetworkTripartition<T>> partitions1, List<NetworkTripartition<T>> partitions2)
    {
        // Count the number tripartitions in the models.
        int ncount = 0, pcount = 0;

        for (NetworkTripartition<T> ntp : partitions1) {
            ncount += ntp.getTripartitionNode().getIndeg();
        }
        for (NetworkTripartition<T> ntp : partitions2) {
            pcount += ntp.getTripartitionNode().getIndeg();
        }

        double fn = (ncount == 0) ? 0.0 : (double) computeTripartitionDiff(partitions1, partitions2) / ncount;
        double fp = (pcount == 0) ? 0.0 : (double) computeTripartitionDiff(partitions2, partitions1) / pcount;
        double avg = (fn + fp) / 2;

        return new double[] {fn, fp, avg};
    }



    /**
     * This function computes the number of tripartitions in <code>models</code> that are not in <code>refs</code>.
     */
    private static <T> int computeTripartitionDiff(List<NetworkTripartition<T>> models, List<NetworkTripartition<T>> refs)
    {
        int diff = 0;
        for (NetworkTripartition<T> ntp : models) {
            if (!refs.contains(ntp)) {
                diff++;
            }
        }

        return diff;
    }



    /**
     * This function computes the normalized tree-based distance between two networks
     *
     * @param net1, net2: Two networks to be compared.
     *
     * @return [fasle-negative, false-positive, average]
     */
    public static <T> double[] computeNormalizedTreeDistance(Network<T> net1, Network<T> net2)
    {
        // Get the trees in the two networks.
        List<Tree> trees1 = new LinkedList<Tree>();	// Trees in net1.
        List<Tree> trees2 = new LinkedList<Tree>();	// Trees in net2.

        for (NetworkTree<T> nt : getTrees(net1)) {
            trees1.add(nt.makeTree());
        }
        for (NetworkTree<T> nt : getTrees(net2)) {
            trees2.add(nt.makeTree());
        }

        return computeNormalizedTreeDistance(trees1, trees2);
    }



    /**
     * This function computes the normalized tree-based distance between two networks given two lists of contained trees
     *
     * @param trees1, trees2: Two lists of contained trees of the two networks under comparison
     *
     * @return [fasle-negative, false-positive, average]
     */
    public static <T> double[] computeNormalizedTreeDistance(List<Tree> trees1, List<Tree> trees2)
    {
        // Initialize a bipartite graph with nodes corresponding to trees in the two networks.
        BipartiteGraph fnBG = new BipartiteGraph(trees1.size(), trees2.size());
        BipartiteGraph fpBG = new BipartiteGraph(trees1.size(), trees2.size());
        BipartiteGraph avgBG = new BipartiteGraph(trees1.size(), trees2.size());

        for (int l = 0; l < trees1.size(); l++) {
            for (int r = 0; r < trees2.size(); r++) {
                SymmetricDifference sd = new SymmetricDifference();
                double fn, fp;
                int ncount, pcount;
                double avg;

                sd.computeDifference(trees1.get(l), trees2.get(r), false);
                fn = sd.getFalseNegativeCount();	// Number of disagreement edges in a left tree.
                fp = sd.getFalsePositiveCount();	// Number of disagreement edges in a right tree.

                ncount = trees1.get(l).getNodeCount() - trees1.get(l).getLeafCount() - 1;	// #Internal edges in ltree.
                pcount = trees2.get(r).getNodeCount() - trees2.get(r).getLeafCount() - 1;	// #Internal edges in rtree.

                fn = (ncount == 0) ? 0.0 : fn / ncount;
                fp = (pcount == 0) ? 0.0 : fp / pcount;
                avg = (fn + fp) / 2.0;

                fnBG.addEdge(l, r, fn);
                fpBG.addEdge(l, r, fp);
                avgBG.addEdge(l, r, avg);
            }
        }

        // Compute the minimum edge covers, which are also the distance between the two networks.
        double fnDist = fnBG.getMinEdgeCoverWeight() / fnBG.getMinEdgeCoverSize();
        double fpDist = fpBG.getMinEdgeCoverWeight() / fpBG.getMinEdgeCoverSize();
        double avgDist = avgBG.getMinEdgeCoverWeight() / avgBG.getMinEdgeCoverSize();

        return new double[] {fnDist, fpDist, avgDist};
    }


    /**
     * This function computes the tree-based distance between two networks.
     *
     * @param net1, net2: The two networks to be compared.
     *
     * @return [fasle-negative, false-positive, average]
     */
    public static <T> double[] computeTreeDistance(Network<T> net1, Network<T> net2)
    {
        // Get the trees in the two networks.
        List<Tree> trees1 = new LinkedList<Tree>();	// Trees in net1.
        List<Tree> trees2 = new LinkedList<Tree>();	// Trees in net2.

        for (NetworkTree<T> nt : getTrees(net1)) {
            trees1.add(nt.makeTree());
        }
        for (NetworkTree<T> nt : getTrees(net2)) {
            trees2.add(nt.makeTree());
        }

        return computeTreeDistance(trees1, trees2);
    }


    /**
     * This function computes the tree-based distance between two networks given the two lists of contained trees
     * as two lists of trees.
     *
     * @param trees1 tree2 lists of trees induced by two networks.
     *
     * @return [fasle-negative, false-positive, average]
     */
    public static <T> double[] computeTreeDistance(List<Tree> trees1, List<Tree> trees2)
    {
        // Initialize a bipartite graph with nodes corresponding to trees in the two networks.
        BipartiteGraph fnBG = new BipartiteGraph(trees1.size(), trees2.size());
        BipartiteGraph fpBG = new BipartiteGraph(trees1.size(), trees2.size());
        BipartiteGraph avgBG = new BipartiteGraph(trees1.size(), trees2.size());

        for (int l = 0; l < trees1.size(); l++) {
            for (int r = 0; r < trees2.size(); r++) {
                SymmetricDifference sd = new SymmetricDifference();
                double fn, fp;
                int ncount, pcount;
                double avg;

                sd.computeDifference(trees1.get(l), trees2.get(r), false);
                fn = sd.getFalseNegativeCount();	// Number of disagreement edges in a left tree.
                fp = sd.getFalsePositiveCount();	// Number of disagreement edges in a right tree.

                ncount = trees1.get(l).getNodeCount() - trees1.get(l).getLeafCount() - 1;	// #Internal edges in ltree.
                pcount = trees2.get(r).getNodeCount() - trees2.get(r).getLeafCount() - 1;	// #Internal edges in rtree.

                fn = (ncount == 0) ? 0.0 : fn / ncount;
                fp = (pcount == 0) ? 0.0 : fp / pcount;
                avg = (fn + fp) / 2.0;

                fnBG.addEdge(l, r, fn);
                fpBG.addEdge(l, r, fp);
                avgBG.addEdge(l, r, avg);
            }
        }

        // Compute the minimum edge covers, which are also the distance between the two networks.
        double fnDist = fnBG.getMinEdgeCoverWeight();
        double fpDist = fpBG.getMinEdgeCoverWeight();
        double avgDist = avgBG.getMinEdgeCoverWeight();

        return new double[] {fnDist, fpDist, avgDist};
    }


    /**
     * This function computes the parsimony for a network N and a list of sequences.
     *
     * @param net: The network to compute the parsimony.
     * @param seq: Raw sequence input string.
     * @param bsize: The block size.
     *
     * @return: The parsimony of <code>net</code> and <code>seqs</code>. It returns -1 if
     * there's any errors during the computation of the parsimony.
     */
    public static <T> int computeParsimony(Network<T> net, String seq, int bsize)
    {
        SequenceAlignment sa = new SequenceAlignment();

        // Get the sequences from the input stream.
        try {
            sa.readSequences(seq);
        }
        catch (SequenceException e) {
            System.err.println("Cannot parse the sequences.");
            System.err.println(e.getMessage());
            return -1;
        }

        // Ready to compute the parsimony.
        return computeParsimony(net, sa, bsize);
    }


    /**
     * This function computes the parsimony for a network N and and a list of sequences.
     *
     * @param net, sa: The network and sequence alignment to compute the parsimony.
     * @param bsize: The block size.
     *
     * @return: The parsimony of <code>net</code> and <code>seqs</code>. It returns -1 if
     * if bsize <= 0, or there's any errors when calling the function computeParsimony of the
     * class ParsimonyCalculator from the phylonet.fitchpars package.
     */
    public static <T> int computeParsimony(Network<T> net, SequenceAlignment sa, int bsize)
    {
        if (bsize <= 0) {
            System.err.println("The block size must be positive");
            return -1;
        }
        else {
            // Get the set of trees induced from the network.
            List<Tree> trees = new LinkedList<Tree>();

            for (NetworkTree<T> nt : getTrees(net)) {
                trees.add(nt.makeTree());
            }

            // Sum the minimum score (over all trees) for all blocks.
            int sum = 0, start = 0;
            ParsimonyCalculator pc = new ParsimonyCalculator();

            while (start < sa.getSequenceLength()) {
                // Get the minimum score for each block.
                int min = Integer.MAX_VALUE;
                int stop = (start + bsize <= sa.getSequenceLength()) ? (start + bsize) : sa.getSequenceLength();

                for (Tree tr : trees) {
                    int score = pc.computeParsimony((MutableTree) tr, sa.getTaxa(), sa.getBlock(start, stop));

                    if (min > score) {
                        min = score;
                    }
                }

                // Take the sum of minimum scores.
                sum += min;

                // Go on to the next block.
                start = stop;
            }

            return sum;
        }
    }


    /**
     * This function checks if two networks have the same topology
     *
     * @param net1, net2: The two networks under comparison
     */
    public static <T> boolean hasTheSameTopology(Network<T> net1, Network<T> net2)
    {
        NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
        return metric.computeDistanceBetweenTwoNetworks(net1,net2)==0;
    }

    public static <T> double computeDistanceBetweenTwoNetworks(Network<T> net1, Network<T> net2)
    {
        NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
        return metric.computeDistanceBetweenTwoNetworks(net1,net2);
    }

    public static <T> String getTopologyString(Network<T> net) {
        Network network = net.clone();
        for(Object nodeObject : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, NetNode.NO_SUPPORT);
                node.setParentDistance(parent, NetNode.NO_DISTANCE);
                node.setParentProbability(parent, NetNode.NO_PROBABILITY);
            }
        }
        return network.toString();
    }

    public static <T> String getCoalUnitString(Network<T> net) {
        Network network = net.clone();
        double popsize = network.getRoot().getRootPopSize();
        for(Object nodeObject : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentDistance(parent, node.getParentDistance(parent) / (popsize / 2));
                node.setParentSupport(parent, NetNode.NO_SUPPORT);
            }
        }
        return network.toString();
    }

    public static <T> String getDendroscopeCompatibleString(Network<T> net) {
        Network network = net.clone();
        for(Object nodeObject : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, NetNode.NO_SUPPORT);
                node.setParentProbability(parent, NetNode.NO_PROBABILITY);
            }
        }
        return network.toString();
    }

    public static <T> String getFullString(Network<T> net) {
        return (!Double.isNaN(net.getRoot().getRootPopSize()) ? ("[" + net.getRoot().getRootPopSize() + "]") : "") + net.toString();
    }


    /**
     * This function adds <code>numReticulations</code> random reticulations to a given network
     * Note that the function now adds reticulations without dealing with branch lengths or inheritance probabilities
     *
     * @param network:  the given network
     * @param numReticulations:     the number of reticulations to add
     */
    public static <T> void addRandomReticulationEdge(Network<T> network, int numReticulations){
        for(int i=0; i<numReticulations; i++){
            addRandomReticulationEdge(network);
        }
    }


    /**
     * This function adds one random reticulation to a given network
     * Note that the function now adds reticulations without dealing with branch lengths or inheritance probabilities
     */
    private static <T> void addRandomReticulationEdge(Network<T> network){
        List<Tuple<NetNode<T>,NetNode<T>>> edgeList = new ArrayList<Tuple<NetNode<T>,NetNode<T>>>();
        Map<NetNode<T>, Integer> node2id = new HashMap<NetNode<T>, Integer>();
        int index = 0;
        for(NetNode<T> node: postTraversal(network)){
            for(NetNode<T> child: node.getChildren()){
                edgeList.add(new Tuple<NetNode<T>, NetNode<T>>(node, child));
            }
            node2id.put(node, index++);
        }
        int size = node2id.size();
        boolean[][] matrix = new boolean[size][size];
        for(NetNode<T> node: postTraversal(network)){
            int nodeID = node2id.get(node);
            matrix[nodeID][nodeID] = true;
            for(NetNode<T> child: node.getChildren()){
                int childID = node2id.get(child);
                matrix[nodeID][childID] = true;
                for(int i=0; i<size; i++){
                    if(matrix[childID][i]){
                        matrix[nodeID][i] = true;
                    }
                }
            }
        }

        int numEdges = edgeList.size();
        int sourceEdgeId = (int)(Math.random() * numEdges);
        Tuple<NetNode<T>,NetNode<T>> sourceEdge = edgeList.get(sourceEdgeId);
        int sourceEdgeChildID = node2id.get(sourceEdge.Item2);
        Tuple<NetNode<T>,NetNode<T>> destinationEdge;
        do{
            destinationEdge = edgeList.get((int)(Math.random() * numEdges));
        }while(matrix[node2id.get(destinationEdge.Item2)][sourceEdgeChildID]);
        NetNode<T> insertedSourceNode = new BniNetNode<T>();
        insertedSourceNode.adoptChild(sourceEdge.Item2, NetNode.NO_DISTANCE);
        sourceEdge.Item1.removeChild(sourceEdge.Item2);
        sourceEdge.Item1.adoptChild(insertedSourceNode, NetNode.NO_DISTANCE);
        NetNode<T> insertedDestinationNode = new BniNetNode<T>();
        insertedDestinationNode.adoptChild(destinationEdge.Item2, NetNode.NO_DISTANCE);
        destinationEdge.Item1.removeChild(destinationEdge.Item2);
        destinationEdge.Item1.adoptChild(insertedDestinationNode, NetNode.NO_DISTANCE);
        insertedSourceNode.adoptChild(insertedDestinationNode, NetNode.NO_DISTANCE);
    }



    /**
     * This function reads a network from its string representation
     */
    public static <T> Network<T> readNetwork(String networkExp){
        try{
            RichNewickReaderAST reader = new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
            reader.setHybridSumTolerance(BigDecimal.valueOf(0.0001));
            NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            RichNewickReadResult<edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks> readResult = reader.read(new ByteArrayInputStream(networkExp.getBytes()));
            return  transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return null;
    }

    /**
     * This function reads a network with root pop from its string representation
     */
    public static <T> Network<T> readNetworkWithRootPop(String networkExp){
        if(networkExp.startsWith("[")) {
            double popSize = Double.parseDouble(networkExp.substring(1, networkExp.indexOf("]")));
            String s = networkExp.substring(networkExp.indexOf("]") + 1);
            Network<T> network = readNetwork(s);
            network.getRoot().setRootPopSize(popSize);
            return network;
        } else {
            return readNetwork(networkExp);
        }
    }

    public static List<String> getTaxaNamesUnderReticulation(Network net) {
        int count = 0;
        List<String> results = new ArrayList<>();
        for(Object leafObject : net.getLeaves()) {
            NetNode leaf = (NetNode) leafObject;
            NetNode node = leaf;
            while(!node.isRoot()) {
                if(node.isNetworkNode()) {
                    results.add(leaf.getName());
                    break;
                }
                node = (NetNode) node.getParents().iterator().next();
            }
        }
        return results;
    }


    /**
     * This function removes all branch lengths and inheritance probabilities from a given network
     */
    public static <T> void removeAllParameters (Network<T> network){
        for(NetNode<T> parent: network.bfs()){
            for(NetNode<T> child: parent.getChildren()){
                child.setParentDistance(parent, NetNode.NO_DISTANCE);
                child.setParentProbability(parent, NetNode.NO_PROBABILITY);
                child.setParentSupport(parent, NetNode.NO_SUPPORT);
            }
        }
    }



    /**
     * This function computes the set of lowest articulation nodes in a network who has at least one child node that is not articulation node
     * @param net: The network
     * @return: The set of articulation nodes
     */
    public static <T> Set<NetNode<T>> getLowestArticulationNodes(Network<T> net){
        List<NetNode<T>> allArticulationNodes = new ArrayList<NetNode<T>>();
        Set<NetNode<T>> articulationNodesToReturn = new HashSet<NetNode<T>>();
        for(NetNode<T> node: Networks.postTraversal(net)){
            if(node.isLeaf()){
                allArticulationNodes.add(node);
            }

            else if(node.isRoot()){
                boolean fArticulate = true;
                for(NetNode child: node.getChildren()){
                    if(!allArticulationNodes.contains(child)){
                        fArticulate = false;
                        break;
                    }
                }
                if(!fArticulate){
                    articulationNodesToReturn.add(node);
                }

            }
            else if(node.isTreeNode()){
                boolean ftotal = true;
                for(NetNode<T> child: node.getChildren()){
                    if(!allArticulationNodes.contains(child)){
                        ftotal = false;
                        break;
                    }
                }
                if(ftotal){
                    allArticulationNodes.add(node);
                }else{
                    NetNode<T> parent = node.getParents().iterator().next();
                    double distance = node.getParentDistance(parent);
                    double support = node.getParentSupport(parent);
                    parent.removeChild(node);
                    boolean disconnect = isDisconnectedNetwork(net, parent);
                    parent.adoptChild(node, distance);
                    node.setParentSupport(parent, support);
                    if (disconnect) {
                        articulationNodesToReturn.add(node);
                        allArticulationNodes.add(node);
                    }
                }

            }
        }
        return articulationNodesToReturn;
    }


    /**
     * This function computes the set of all articulation nodes in a network
     * @param net: The network
     * @return: The set of articulation nodes
     */
    public static <T> Set<NetNode<T>> getAllArticulationNodes(Network<T> net){
        Set<NetNode<T>> articulationNodes = new HashSet<>();
        for(NetNode<T> node: Networks.postTraversal(net)){
            if(node.isLeaf() || node.isRoot()){
                articulationNodes.add(node);
            }
            else if(node.isTreeNode()){
                boolean childrenArticulate = true;
                for(NetNode<T> child: node.getChildren()){
                    if(!articulationNodes.contains(child)){
                        childrenArticulate = false;
                        break;
                    }
                }
                if(childrenArticulate){
                    articulationNodes.add(node);
                }else{
                    NetNode<T> parent = node.getParents().iterator().next();
                    double distance = node.getParentDistance(parent);
                    double support = node.getParentSupport(parent);
                    parent.removeChild(node);
                    boolean disconnect = isDisconnectedNetwork(net, parent);
                    parent.adoptChild(node, distance);
                    node.setParentSupport(parent, support);
                    if (disconnect) {
                        articulationNodes.add(node);
                    }
                }

            }
        }
        return articulationNodes;
    }



    /**
     * This function is to help compute articulate nodes
     */
    public static <T> boolean isDisconnectedNetwork(Network<T> net, NetNode<T> ignoreNode){
        Set<NetNode<T>> visited = new HashSet<NetNode<T>>();
        Set<NetNode<T>> seen = new HashSet<NetNode<T>>();
        for(NetNode<T> node: Networks.postTraversal(net)){
            if(node.getIndeg()==1 && node.getOutdeg()==1 && node!=ignoreNode){
                return false;
            }
            visited.add(node);
            for(NetNode<T> parent: node.getParents()){
                seen.add(parent);
            }
            for(NetNode<T> child: node.getChildren()){
                seen.add(child);
            }
        }

        return visited.size()==seen.size();
    }


    /**
     * This function returns nodes in a network in post traversal order
     */
    public static <T> List<NetNode<T>> postTraversal(Network<T> net){
        Stack<NetNode<T>> stack = new Stack<NetNode<T>>();
        List<NetNode<T>> searchedNodes = new ArrayList<NetNode<T>>();
        stack.push(net.getRoot());
        Map<NetNode<T>, Integer> node2index = new HashMap<NetNode<T>, Integer>();
        node2index.put(net.getRoot(), 0);

        while(!stack.isEmpty()){
            NetNode<T> topNode = stack.peek();
            int index = node2index.get(topNode);
            if(index == topNode.getOutdeg()){
                searchedNodes.add(stack.pop());
            }
            else{
                Iterator<NetNode<T>> it = topNode.getChildren().iterator();
                for(int i=0; i<index; i++){
                    it.next();
                }
                NetNode<T> child = it.next();
                if(searchedNodes.contains(child)){
                    node2index.put(topNode, index + 1);
                }
                else{
                    stack.push(child);
                    node2index.put(child, 0);
                }
            }
        }
        return searchedNodes;
    }


    /**
     * This function checks if a network contains cycle
     */
    public static <T> boolean hasCycle(Network<T> net) {
        NetNode<T> root = net.getRoot();
        Set<NetNode<T>> visited = new HashSet<NetNode<T>>();
        Set<NetNode<T>> stacked = new HashSet<NetNode<T>>();
        if(dfs(root, visited, stacked)) return true;
        return false;
    }

    public static <T> boolean hasParallelEdges(Network<T> net) {
        for(NetNode<T> node : net.bfs()) {
            if(node.isNetworkNode()) {
                Iterator<NetNode<T>> it = node.getParents().iterator();
                NetNode<T> parent1 = it.next();
                NetNode<T> parent2 = it.next();
                if (parent1 == parent2) return true;
            }
        }
        return false;
    }


    /**
     * This function returns nodes in a network in dfs order
     */
    private static <T> boolean dfs(NetNode<T> node, Set<NetNode<T>> visited,
                                   Set<NetNode<T>> stacked) {
        if(stacked.contains(node)) return true;
        if(visited.contains(node)) return false;

        visited.add(node);
        stacked.add(node);
        for(NetNode<T> n : node.getChildren()) {
            if(dfs(n, visited, stacked)) return true;
        }
        stacked.remove(node);
        return false;
    }


    /**
     * This function computes the maximal load of a network
     */
    public static <T> double computeMaxLoad(Network<T> network, Map<String, List<String>> species2alleles){
        NetworkLoad networkLoad = new NetworkLoad();
        return networkLoad.computeMaxLoad(network, species2alleles);
    }


    /**
     * Map nodes in two networks
     * @param net1   network 1
     * @param net2   network 2
     * @param <T>
     * @return       map result
     */
    public static <T> Map<NetNode<T>, NetNode<T>> mapTwoNetworks(Network<T> net1, Network<T> net2) {
        Map<NetNode<T>, NetNode<T>> map = new HashMap<NetNode<T>, NetNode<T>>();
        if(mapTwoNodes(net1.getRoot(), net2.getRoot(), map)) return map;
        else return null;
    }

    /**
     * Recursively map two nodes
     * @param node1   node 1
     * @param node2   node 2
     * @param map     final map
     * @param <T>
     * @return        map successfully or not
     */
    private static <T> boolean mapTwoNodes(NetNode<T> node1, NetNode<T> node2,
                                          Map<NetNode<T>, NetNode<T>> map) {
        if(node1.isLeaf() && node2.isLeaf()) {
            if(node1.getName().equals(node2.getName())) {
                map.put(node1, node2);
                return true;
            }
            return false;
        } else if(node1.isLeaf() || node2.isLeaf()) {
            return false;
        }
        // non-leave nodes
        List<NetNode<T>> children1list = IterableHelp.toList(node1.getChildren());
        Set<NetNode<T>> children1 = new HashSet<NetNode<T>>(children1list);
        List<NetNode<T>> children2list = IterableHelp.toList(node2.getChildren());
        Set<NetNode<T>> children2 = new HashSet<NetNode<T>>(children2list);
        if(children1.size() != children2.size()) return false;

        // perform mapping for each children
        for(NetNode<T> c1 : children1) {
            // if already in map
            if(map.containsKey(c1)) {
                if(children2.contains(map.get(c1))) {
                    children2.remove(map.get(c1));
                    continue;
                } else return false;
            }
            // if not, try to map to one of the nodes in children2
            NetNode<T> mapNode2 = null;
            for(NetNode<T> c2 : children2) {
                if(mapTwoNodes(c1, c2, map)) {
                    mapNode2 = c2;
                    break;
                }
            }
            if(mapNode2 == null) return false; // false to map c1 to any node in children2
            children2.remove(mapNode2);
        }
        map.put(node1, node2);
        return true;
    }

    /**
     * Gets the edges in a network
     * @param net   the network
     * @param <T>   data type
     * @return  all the edges
     */
    public static <T> List<Tuple<NetNode<T>, NetNode<T>>> getAllEdges(Network<T> net) {
        List<Tuple<NetNode<T>, NetNode<T>>> edges = new ArrayList<>();
        for(NetNode<T> node : net.dfs()) {
            for(NetNode<T> par : node.getParents()) {
                edges.add(new Tuple<>(node, par));
            }
        }
        return edges;
    }


    /**
     * Get the other child of a network tree node given a child
     * @param parent    the network tree node
     * @param child     one of the child
     * @param <T>
     * @return  the other child
     */
    public static <T> NetNode<T> getOtherChild(NetNode<T> parent, NetNode<T> child) {
        boolean found = false;
        NetNode<T> otherChild = null;
        for(NetNode<T> node : parent.getChildren()) {
            if(node == child) {
                found = true;
            } else {
                otherChild = node;
            }
        }
        if(parent.getChildCount() != 2 || !found || otherChild == null) {
            throw new RuntimeException("Invalid node with wrong in/out-degree " + parent.getName());
        }
        return otherChild;
    }


    /**
     * Get the other parent of a network node
     * @param child     the network node
     * @param parent    one of the parents
     * @param <T>
     * @return  the other parent
     */
    public static <T> NetNode<T> getOtherParent(NetNode<T> child, NetNode<T> parent) {
        boolean found = false;
        NetNode<T> otherParent = null;
        for(NetNode<T> node : child.getParents()) {
            if(node == parent) {
                found = true;
            } else {
                otherParent = node;
            }
        }
        if(child.getParentCount() != 2 || !found || otherParent == null) {
            throw new RuntimeException("Invalid node with wrong in/out-degree " + child.getName());
        }
        return otherParent;
    }


    /**
     * Get all internal nodes of a network
     * @param net   the network
     * @param <T>
     * @return  all internal nodes
     */
    public static <T> List<NetNode<T>> getInternalNodes(Network<T> net) {
        List<NetNode<T>> res = new ArrayList<>();
        Queue<NetNode<T>> queue = new LinkedList<>();
        queue.add(net.getRoot());
        Set<NetNode<T>> visited = new HashSet<>();
        while(!queue.isEmpty()) {
            NetNode<T> node = queue.poll();
            if(visited.contains(node)) continue;
            res.add(node);
            visited.add(node);
            for(NetNode<T> child: node.getChildren()) {
                if(child.isLeaf()) continue;
                queue.add(child);
            }
        }
        return res;
    }


    /**
     * Get all internal tree nodes
     * @param net   the network
     * @param <T>
     * @return  all internal tree nodes
     */
    public static <T> List<NetNode<T>> getInternalTreeNodes(Network<T> net) {
        List<NetNode<T>> res = new ArrayList<>();
        Queue<NetNode<T>> queue = new LinkedList<>();
        queue.add(net.getRoot());
        Set<NetNode<T>> visited = new HashSet<>();
        while(!queue.isEmpty()) {
            NetNode<T> node = queue.poll();
            if(visited.contains(node)) continue;
            visited.add(node);
            if(node.isTreeNode()) res.add(node);
            for(NetNode<T> child: node.getChildren()) {
                if(child.isLeaf()) continue;
                queue.add(child);
            }
        }
        return res;
    }

    public static <T> void computeBiconComponetEdges(List<List<Tuple<NetNode,NetNode>>> results,
                                                     Stack<Tuple<NetNode,NetNode>> S,
                                                     Map<NetNode, Integer> num,
                                                     Map<NetNode, Integer> lowpt,
                                                     Integer index,
                                                     NetNode<T> v,
                                                     NetNode<T> u) {
        index = index + 1;
        num.put(v, index);
        lowpt.put(v, index);
        List<NetNode<T>> adj = new ArrayList<>();
        for(NetNode<T> node : v.getChildren())
            adj.add(node);
        for(NetNode<T> node : v.getParents())
            adj.add(node);
        for(NetNode<T> w : adj) {
            Tuple<NetNode,NetNode> e = new Tuple(v, w);
            if(num.get(w) == 0) {
                S.push(e);
                computeBiconComponetEdges(results, S, num, lowpt, index, w, v);
                lowpt.put(v, Math.min(lowpt.get(v), lowpt.get(w)));
                if(lowpt.get(w) >= num.get(v)) {
                    List<Tuple<NetNode,NetNode>> result = new ArrayList<>();
                    Tuple<NetNode,NetNode> curE;
                    do {
                        curE = S.pop();
                        result.add(curE);
                    } while(!curE.equals(e));
                    results.add(result);
                }
            } else if(num.get(w) < num.get(v) && w != u){
                S.push(e);
                lowpt.put(v, Math.min(lowpt.get(v), num.get(w)));
            }
        }
    }

    public static <T> List<List<Tuple<NetNode,NetNode>>> computeBiconComponets(Network<T> network) {
        Integer index = 0;
        Stack<Tuple<NetNode,NetNode>> S = new Stack<>();
        Map<NetNode, Integer> num = new HashMap<>();
        Map<NetNode, Integer> lowpt = new HashMap<>();
        for(NetNode node : network.dfs()) {
            num.put(node, 0);
            lowpt.put(node, 0);
        }
        List<List<Tuple<NetNode,NetNode>>> results = new ArrayList<>();
        for(NetNode node : network.dfs()) {
            if(num.get(node) == 0)
                computeBiconComponetEdges(results, S, num, lowpt, index, node, null);
        }

        return results;
    }

    public static <T> int computeLevel(Network<T> network) {
        List<List<Tuple<NetNode,NetNode>>> bicons = computeBiconComponets(network);
        int num = 0;
        for(List<Tuple<NetNode,NetNode>> edges : bicons) {
            Set<NetNode> reticulateNodes = new HashSet<>();
            for(Tuple<NetNode,NetNode> e : edges) {
                if(e.Item1.isNetworkNode())
                    reticulateNodes.add(e.Item1);
                if(e.Item2.isNetworkNode())
                    reticulateNodes.add(e.Item2);
            }
            num = Math.max(num, reticulateNodes.size());
        }
        return num;

    }

    public static <T> boolean isNestedNetwork(Network<T> network) {
        for(NetNode node : network.getNetworkNodes()) {
            for(Object obj : node.getParents()) {
                NetNode parent = (NetNode) obj;

                while(!parent.isRoot()) {
                    if(parent.isNetworkNode())
                        return false;
                    parent = (NetNode) parent.getParents().iterator().next();
                }
            }
        }
        return true;
    }

    public static <T> boolean isGalledNetwork(Network<T> network) {
        for(NetNode node : network.getNetworkNodes()) {
            Iterator<NetNode> parent = node.getParents().iterator();
            NetNode parent1 = parent.next();
            NetNode parent2 = parent.next();

            Set<NetNode> visited = new HashSet<>();
            while(!parent1.isNetworkNode()) {
                visited.add(parent1);
                if(parent1.isRoot()) break;
                parent1 = (NetNode) parent1.getParents().iterator().next();
            }
            boolean flag = false;
            while(!parent2.isNetworkNode()) {
                if(visited.contains(parent2)) {
                    flag = true;
                    break;
                }
                if(parent2.isRoot()) break;
                parent2 = (NetNode) parent2.getParents().iterator().next();
            }
            if(!flag) {
                return false;
            }
        }
        return true;
    }


    public static <T> void scaleNetwork(Network<T> net, double scale) {
        for(NetNode<T> node : postTraversal(net)) {
            for(NetNode<T> child : node.getChildren()) {
                child.setParentDistance(node, child.getParentDistance(node) * scale);
            }
        }
    }


    // Data members
    public static final String NAME_PREFIX = "I";	// Name prefix for interior nodes.

    /**
     * @Description: remove names of all the internal nodes in network net
     * @Param: net: The network to remove names
     * @Author: Zhen Cao
     * @Date: 2019-08-15
     */
    public static <T> void removeInternalNodeNames(Network<T> net){
        for(NetNode<T> node : postTraversal(net)) {
            if (!node.isLeaf()){
                node.setName("");
            }
        }
    }

    /**
     * @return <code>true</code> if these networks have identical leafsets
     * by name.
     */
    public static final <T> boolean leafSetsAgree(Network<T> network1, Network<T> network2) {
        // They definitely can't agree if they are different sizes
        if (network1.getLeafCount() != network2.getLeafCount()) {
            return false;
        }

        Set<String> network1Leaves = new HashSet<>();
        Set<String> network2Leaves = new HashSet<>();

        network1.getLeaves().forEach(leaf -> network1Leaves.add(leaf.getName()));
        network2.getLeaves().forEach(leaf -> network2Leaves.add(leaf.getName()));

        return network1Leaves.containsAll(network2Leaves) && network2Leaves.containsAll(network1Leaves);
    }
}