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
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.bipartitematching.HungarianBipartiteMatcher;
import edu.rice.cs.bioinfo.programs.phylonet.algos.fitchpars.ParsimonyCalculator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTreeEnumerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTripartition;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.sequence.model.SequenceAlignment;
import edu.rice.cs.bioinfo.programs.phylonet.structs.sequence.model.SequenceException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.ByteArrayInputStream;
import java.io.StringWriter;
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

/*
	public static <T> void removeBinaryNodes(Network<T> net)
	{
		// Find all binary nodes.
		List<NetNode<T>> binary_nodes = new LinkedList<NetNode<T>>();
		for (NetNode<T> node : net.bfs()) {
            //System.out.println(node.getName() + " " + node.getIndeg() + " " + node.getOutdeg());
			if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
				binary_nodes.add(node);
			}
		}

		// Remove them.
		for (NetNode<T> node : binary_nodes) {
			NetNode<T> parent = node.getParents().iterator().next();	// Node's only parent.
			NetNode<T> child = node.getChildren().iterator().next();	// Node's only child.
			double distance = node.getParentDistance(parent) + child.getParentDistance(node);
			double inheritanceProb = child.getParentProbability(node);
            //System.out.println("removing " + node.getName());
			parent.removeChild(node);
            //System.out.println("removing " + child.getName());
			node.removeChild(child);
			parent.adoptChild(child, distance);
            child.setParentProbability(parent, inheritanceProb);
		}
	}
    */


    /**
     * This function removes all binary nodes from the network.
     */
    public static <T> void removeBinaryNodes(Network<T> net) {
        boolean update;
        do{
            update = false;
            for (NetNode<T> node : net.bfs()) {
                //System.out.println(node.getName() + " " + node.getIndeg() + " " + node.getOutdeg());
                if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
                    NetNode<T> parent = node.getParents().iterator().next();    // Node's only parent.
                    NetNode<T> child = node.getChildren().iterator().next();    // Node's only child.
                    double distance = node.getParentDistance(parent) + child.getParentDistance(node);
                    double inheritanceProb = child.getParentProbability(node);
                    //System.out.println("removing " + node.getName());
                    parent.removeChild(node);
                    //System.out.println("removing " + child.getName());
                    node.removeChild(child);
                    parent.adoptChild(child, distance);
                    child.setParentProbability(parent, inheritanceProb);
                    update = true;
                }
				if(node.isRoot() && node.getOutdeg()==1){
					NetNode newRoot = node.getChildren().iterator().next();
					node.removeChild(newRoot);
					net.resetRoot(newRoot);
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
	 * This function compute the cluster-based distance between two networks.
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


	public static <T> double[] computeSoftwiredClusterDistance(List<NetworkCluster<T>> clusters1, List<NetworkCluster<T>> clusters2)
	{
		double fn = (clusters1.size() == 0) ? 0.0 : (double) computeSoftwiredClusterDiff(clusters1, clusters2) / clusters1.size();
		double fp = (clusters2.size() == 0) ? 0.0 : (double) computeSoftwiredClusterDiff(clusters2, clusters1) / clusters2.size();
		double avg = (fn + fp) / 2;

		return new double[] {fn, fp, avg};
	}


    /**
     * This function compute the cluster-based distance between two networks.
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


    public static <T> double[] computeHardwiredClusterDistance(List<Set<String>> clusters1, List<Set<String>> clusters2)
    {
        double fn = (clusters1.size() == 0) ? 0.0 : (double) computeHardwiredClusterDiff(clusters1, clusters2) / clusters1.size();
        double fp = (clusters2.size() == 0) ? 0.0 : (double) computeHardwiredClusterDiff(clusters2, clusters1) / clusters2.size();
        double avg = (fn + fp) / 2;

        return new double[] {fn, fp, avg};
    }

    /*
    private static Set<String> convertHardwiredClustersToStringRepresentatives(List<Set<String>> clusters){
        Set<String> representatives = new HashSet<String>();
        for(Set<String> cluster: clusters){
            Arrays.so
        }
        return representatives;
    }
    */


    /**
	 * This function computes the number of clusters in <code>models</code> that are not in <code>refs</code>.
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
     * This function computes the number of clusters in <code>models</code> that are not in <code>refs</code>.
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
				fp = sd.getFalseNegativeCount();	// Number of disagreement edges in a right tree.

				ncount = trees1.get(l).getNodeCount() - trees1.get(l).getLeafCount() - 1;	// #Internal edges in ltree.
				pcount = trees2.get(r).getNodeCount() - trees2.get(r).getLeafCount() - 1;	// #Internal edges in rtree.

				fn = (ncount == 0) ? 0.0 : fn / ncount;
				fp = (pcount == 0) ? 0.0 : fp / pcount;
				avg = (fn + fp) / 2.0;

				fnBG.addEdge(l, r, fn);
				fpBG.addEdge(l, r, avg);
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
	 * This function computes the tree-based distance between two networks represented
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
				fp = sd.getFalseNegativeCount();	// Number of disagreement edges in a right tree.

				ncount = trees1.get(l).getNodeCount() - trees1.get(l).getLeafCount() - 1;	// #Internal edges in ltree.
				pcount = trees2.get(r).getNodeCount() - trees2.get(r).getLeafCount() - 1;	// #Internal edges in rtree.

				fn = (ncount == 0) ? 0.0 : fn / ncount;
				fp = (pcount == 0) ? 0.0 : fp / pcount;
				avg = (fn + fp) / 2.0;
                
				fnBG.addEdge(l, r, fn);
				fpBG.addEdge(l, r, avg);
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
	 * This function computes the parsimony for a network N and and a list of sequences.
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


    public static <T> boolean hasTheSameTopology(Network<T> net1, Network<T> net2)
    {
        NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
        return metric.computeDistanceBetweenTwoNetworks(net1,net2)==0;
    }

    //With respect to topology only
    public static <T> void addRandomReticulationEdge(Network<T> network, int numReticulations){
        for(int i=0; i<numReticulations; i++){
            addRandomReticulationEdge(network);
        }
    }


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


    public static <T> void removeAllParameters (Network<T> network){
        for(NetNode<T> parent: network.bfs()){
            for(NetNode<T> child: parent.getChildren()){
                child.setParentDistance(parent, NetNode.NO_DISTANCE);
                child.setParentProbability(parent, NetNode.NO_PROBABILITY);
                child.setParentSupport(parent, NetNode.NO_SUPPORT);
            }
        }
    }



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
        //System.out.println("Source Edge:(" + sourceEdge.Item1.getName()+ "," +sourceEdge.Item2.getName()+ ")");
        int sourceEdgeChildID = node2id.get(sourceEdge.Item2);
        Tuple<NetNode<T>,NetNode<T>> destinationEdge;
        do{
            destinationEdge = edgeList.get((int)(Math.random() * numEdges));
            //System.out.println("Destination Edge:(" + destinationEdge.Item1.getName()+ "," +destinationEdge.Item2.getName()+ ")");
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
     * This function computes the set of lowest articulation nodes in a networks whose has at
     * least one child node that is not articulation node
     * @param net, sa: The network
     * * @return: The set of articulation nodes
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
                    parent.removeChild(node);
                    boolean disconnect = isDisconnectedNetwork(net, parent);
                    parent.adoptChild(node, distance);
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
     * This function computes the set of articulation nodes in a networks whose has at
     * least one child node that is not articulation node
     * @param net, sa: The network
     * * @return: The set of articulation nodes
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
                    parent.removeChild(node);
                    boolean disconnect = isDisconnectedNetwork(net, parent);
                    parent.adoptChild(node, distance);
                    if (disconnect) {
                        articulationNodes.add(node);
                    }
                }

            }
        }
        return articulationNodes;
    }






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


    public static <T> boolean hasCycle(Network<T> net) {
        NetNode<T> root = net.getRoot();
        Set<NetNode<T>> visited = new HashSet<NetNode<T>>();
        Set<NetNode<T>> stacked = new HashSet<NetNode<T>>();
        if(dfs(root, visited, stacked)) return true;
        return false;
    }

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



    public static <T> double computeMaxLoad(Network<T> network, Map<String, List<String>> species2alleles){
        NetworkLoad networkLoad = new NetworkLoad();
        return networkLoad.computeMaxLoad(network, species2alleles);
    }


	// Data members
	public static final String NAME_PREFIX = "I";	// Name prefix for interior nodes.
}



/**
 * This class provides a simple way to represent the bipartite graph. The class Networks
 * use this class to compute the tree-based distance between two networks.
 */
class BipartiteGraph {
	/**
	 * This constructor initializes a bipartite graph with <code>lsize</code> nodes on
	 * the left and <code>rsize</code> nodes on the right.
	 */
	public BipartiteGraph(int lsize, int rsize)
	{
		assert(lsize > 0 && rsize > 0);

		_lsize = lsize;
		_rsize = rsize;

		_weights = new double[lsize][rsize];
		for (int l = 0; l < lsize; l++) {
			for (int r = 0; r < rsize; r++) {
				_weights[l][r] = NO_EDGE;
			}
		}
	}

	/**
	 * This function adds an edge to the bipartite graph.
	 *
	 * @param l, r: indexes of the left and right endpoints for this edge.
	 * @param w: the weight of this edge. We require it to be positive.
	 */
	public void addEdge(int l, int r, double w)
	{
		assert(0 <= l && l < _lsize);
		assert(0 <= r && r < _rsize);
		assert(w >= 0);

		_weights[l][r] = w;
	}

	/**
	 * This function returns the weight of a minimum edge cover of this bipartite graph.
	 *
	 * @return: The weight of a minimum edge cover.
	 */
	public double getMinEdgeCoverWeight()
	{
		computeMinEdgeCover();

		double weight = _hbm.getMatchingWeight();
		weight = _hbm.getMatching().length * _invert_weight - weight;

		return weight / 2.0;
	}

	/**
	 * Computes and returns the size of a minimum edge cover of this bipartite graph.
	 *
	 * @return: The number of edges in a minimum edge cover.
	 */
	public int getMinEdgeCoverSize()
	{
		computeMinEdgeCover();

		int matching[][] = _hbm.getMatching();
		int size = matching.length;

		for (int i = 0; i < matching.length; i++) {
			int l = matching[i][0];
			int r = matching[i][1];

			if (l >= _lsize && r < _lsize) {
				size--;	// Edges in the copy of the graph is removed from the cover.
			}
		}

		return size;
	}

	/**
	 * This function initializes the bipartite matcher to find the minimum-weight edge
	 * cover. Refer to Chapters 17 and 19 in the book "Combinatorial Optimzation" by A. Schrijver
	 * to see how the Hungarian algorithm is used to find the minimum edge cover.
	 */
	private void initBipartiteMatcher()
	{
		_invert_weight = 2 * getMaxEdgeWeight() + 1.0;

		// Initialize the original and disjoint-copy parts for the new bipartite graph.
		_hbm = new HungarianBipartiteMatcher(_lsize + _rsize, _lsize + _rsize);

		for (int l = 0; l < _lsize; l++) {
			for (int r = 0; r < _rsize; r++) {
				if (_weights[l][r] != NO_EDGE) {
					_hbm.addEdge(l, r + _lsize, _invert_weight - _weights[l][r]);	// Original part.
					_hbm.addEdge(r + _lsize, l, _invert_weight - _weights[l][r]);	// Disjoint-copy part.
				}
			}
		}

		// Add edges for corresponding nodes.
		for (int l = 0; l < _lsize; l++) {
			double min = getMinIncidentEdgeWeight(l, true);

			_hbm.addEdge(l, l, _invert_weight -  2 * min);
		}
		for (int r = 0; r < _rsize; r++) {
			double min = getMinIncidentEdgeWeight(r, false);

			_hbm.addEdge(r + _lsize, r + _lsize, _invert_weight -  2 * min);
		}
	}

	/**
	 * Computes the minimum edge covers for this bipartite graph.
	 */
	private void computeMinEdgeCover()
	{
		initBipartiteMatcher();
		_hbm.findMatching();
	}

	/**
	 * This function returns the minimum weight of all edges incident with a node indexed by
	 * <code>index</index>.
	 *
	 * @param index: The index of the node, either in the left side or the right side.
	 * @param left: <code>true</code> if this node is on the left side; <code>false</code> otherwise.
	 */
	private double getMinIncidentEdgeWeight(int index, boolean left)
	{
		double min = Double.POSITIVE_INFINITY;

		if (left) {	// Find the min weight for a node on the left side.
			for (int r = 0; r < _rsize; r++) {
				if (min > _weights[index][r] && _weights[index][r] != NO_EDGE) {
					min = _weights[index][r];
				}
			}
		}
		else {	// Find the min weight for a node on the right side.
			for (int l = 0; l < _lsize; l++) {
				if (min > _weights[l][index] && _weights[l][index] != NO_EDGE) {
					min = _weights[l][index];
				}
			}
		}

		return min;
	}

	/**
	 * This function returns the index of the edge with minimum weight of all edges incident
	 * with a node indexed by <code>index</index>.
	 *
	 * @param index: The index of the node, either in the left side or the right side.
	 * @param left: <code>true</code> if this node is on the left side; <code>false</code> otherwise.
	 */
	private int getMinIncidentEdgeIndex(int index, boolean left)
	{
		double min = Double.POSITIVE_INFINITY;
		int i = -1;	// Indicating that this node has no incident edges.

		if (left) {	// Find the min weight for a node on the left side.
			for (int r = 0; r < _rsize; r++) {
				if (min > _weights[index][r] && _weights[index][r] != NO_EDGE) {
					min = _weights[index][r];
					i = r;
				}
			}
		}
		else {	// Find the min weight for a node on the right side.
			for (int l = 0; l < _lsize; l++) {
				if (min > _weights[l][index] && _weights[l][index] != NO_EDGE) {
					min = _weights[l][index];
					i = l;
				}
			}
		}

		return i;
	}

	/**
	 * This function returns the maximum weight of all edges in the bipartite graph.
	 */
	private double getMaxEdgeWeight()
	{
		double max = Double.NEGATIVE_INFINITY;

		for (int l = 0; l < _lsize; l++) {
			for (int r = 0; r < _rsize; r++) {
				if (max < _weights[l][r] && _weights[l][r] != NO_EDGE) {
					max = _weights[l][r];
				}
			}
		}

		return max;
	}




	// Data members
	public final static double NO_EDGE = Double.NEGATIVE_INFINITY;

	private int _lsize;		// The number of nodes on the left side.
	private int _rsize;		// The number of nodes on the right side.
	private double _weights[][];			// Holds the edge weights.

	private HungarianBipartiteMatcher _hbm;	// The Hungarian matcher used to find the max. matching.
	private double _invert_weight;			// Hold a "big enough" weight, allowing the matcher to find min. matching.
}
