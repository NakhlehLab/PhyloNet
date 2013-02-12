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

import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.bipartitematching.HungarianBipartiteMatcher;
import edu.rice.cs.bioinfo.programs.phylonet.algos.fitchpars.ParsimonyCalculator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTreeEnumerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTripartition;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.sequence.model.SequenceAlignment;
import edu.rice.cs.bioinfo.programs.phylonet.structs.sequence.model.SequenceException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

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
	public static <T> void removeBinaryNodes(Network<T> net)
	{
		// Find all binary nodes.
		List<NetNode<T>> binary_nodes = new LinkedList<NetNode<T>>();
		for (NetNode<T> node : net.bfs()) {
			if (node.getIndeg() == 1 && node.getOutdeg() == 1) {
				binary_nodes.add(node);
			}
		}

		// Remove them.
		for (NetNode<T> node : binary_nodes) {
			NetNode<T> parent = node.getParents().iterator().next();	// Node's only parent.
			NetNode<T> child = node.getChildren().iterator().next();	// Node's only child.
			double distance = node.getParentDistance(parent) + child.getParentDistance(node);
			//double gamma = node.getGamma(0) * parent.getGamma(node);
			parent.removeChild(node);
			node.removeChild(child);
		//	parent.adoptChild(child, distance, gamma);
		}
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
	public static <T> Iterable<NetworkCluster<T>> getClusters(Network<T> net)
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
	public static <T> double[] computeClusterDistance(Network<T> net1, Network<T> net2)
	{
		List<NetworkCluster<T>> clusters1 = new LinkedList<NetworkCluster<T>>();
		List<NetworkCluster<T>> clusters2 = new LinkedList<NetworkCluster<T>>();

		for (NetworkCluster<T> nc : getClusters(net1)) {
			clusters1.add(nc);
		}
		for (NetworkCluster<T> nc : getClusters(net2)) {
			clusters2.add(nc);
		}

		return computeClusterDistance(clusters1, clusters2);
	}

	/**
	 * Compute the cluter-based distance between two networks represented as two lists
	 * of clusters.
	 *
	 * @param Two lists of clusters representing two networks.
	 *
	 * @return [fasle-negative, false-positive, average]
	 */
	public static <T> double[] computeClusterDistance(List<NetworkCluster<T>> clusters1, List<NetworkCluster<T>> clusters2)
	{
		double fn = (clusters1.size() == 0) ? 0.0 : (double) computeClusterDiff(clusters1, clusters2) / clusters1.size();
		double fp = (clusters2.size() == 0) ? 0.0 : (double) computeClusterDiff(clusters2, clusters1) / clusters2.size();
		double avg = (fn + fp) / 2;

		return new double[] {fn, fp, avg};
	}

	/**
	 * This function computes the number of clusters in <code>models</code> that are not in <code>refs</code>.
	 */
	private static <T> int computeClusterDiff(List<NetworkCluster<T>> models, List<NetworkCluster<T>> refs)
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
	 * This function computes the distance between two networks represented as two lists of tripartitions.
	 *
	 * @param Two networks represented as two list of tripartitions.
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

				sd.computeDifference(trees1.get(l), trees2.get(r));
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
	 * @param Two lists of trees induced by two networks.
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

				sd.computeDifference(trees1.get(l), trees2.get(r));
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

	/**
	 * Restrict network 1 and network 2 on the common set of leaves.
	 *
	 * @param net1
	 * @param net2
	 */
	public <T> void pruneNetworks(Network<T> net1, Network<T> net2) {

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
	 * @param weight: the weight of this edge. We require it to be positive.
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
