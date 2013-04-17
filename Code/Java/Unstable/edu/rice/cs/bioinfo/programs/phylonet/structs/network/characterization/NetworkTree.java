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

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/28/11
 * Time: 3:14 PM
 * To change this template use File | Settings | File Templates.
 */

public class NetworkTree<T> {
	/**
	 * This constructor initializes a network tree.
	 *
	 * @param net: The network that this tree belongs to.
	 * @param edges: The set of network edges chosen to build this tree.
	 */
	public NetworkTree(Network<T> net, List<NetworkEdge<T>> edges)
	{
		_net = net;
		_edges = edges;
		_clusters = null;
	}

	/**
	 * This class returns a <code>STITree</code> counterpart for this network tree. The tree
	 * returned will *not* contain redundant nodes (i.e. nodes with indeg = outdeg = 1).
	 *
	 * @return The <code>Tree<code> corresponding to this network tree.
	 */
	public Tree makeTree()
	{
		Queue<NetNode<T>> source = new LinkedList<NetNode<T>>();	// Nodes in the network used to build the tree.
		Queue<TMutableNode> dest = new LinkedList<TMutableNode>();	// Nodes in the tree built.
		LinkedList<TMutableNode> orphans = new LinkedList<TMutableNode>();	// List of tree leaves corresponding to network nodes.
		MutableTree tree = new STITree<T>();	// The tree to be built.

		// Initialize the two queues to start building the tree.
		if (tree.isEmpty()) {
			tree.createRoot();
		}
		source.offer(_net.getRoot());
		dest.offer(tree.getRoot());

		// Iterate over the network and build the tree from the list of edges.
		while (!source.isEmpty()) {
			NetNode<T> parent = source.poll();
			TMutableNode peer = dest.poll();	// Corresponding tree node to parent.

			for (NetNode<T> child : parent.getChildren()) {
				if (child.isTreeNode() || _edges.contains(new NetworkEdge<T>(parent, child))) {
					// Create a tree node corresponding to child under peer.
					TMutableNode copy;
					if (child.getName() == NetNode.NO_NAME) {
						copy = peer.createChild(TNode.NO_NAME);
					}
					else {
						copy = peer.createChild(child.getName());
					}

					// Update the distance and data for this child.
					double distance = child.getParentDistance(parent);
					if (distance == NetNode.NO_DISTANCE) {
						copy.setParentDistance(TNode.NO_DISTANCE);
					}
					else {
						copy.setParentDistance(distance);
					}
					((STINode<T>) copy).setData(child.getData());

					// Continue to iterate over the children of nn and tn.
					source.offer(child);
					dest.offer(copy);
				}
			}

			if (peer.getChildCount() == 1) {
				suppressBinaryNode(peer);
			}
			if (!parent.isLeaf() && peer.isLeaf()) {
				orphans.add(peer);
			}

		}

		// Eliminate bad nodes.
		for (TMutableNode node : orphans) {
			removeOrphan(node);
		}

		return tree;
	}

	/**
	 * This function suppresses tree nodes that has only one child.
	 *
	 * This function assumes that node has only one child. This is checked before
	 * calling this function.
	 */
	private void suppressBinaryNode(TMutableNode node)
	{
		assert(node.getChildCount() == 1);

		if (!node.isRoot()) {	// Interior node
			// Move node's sole child up one level to node's parent.
			TMutableNode parent = node.getParent();
			parent.removeChild(node, true);
		}
		else {	// Root
			TMutableNode child = node.getChildren().iterator().next();	// The sole child.
			child.makeRoot();
			child.removeChild(node, false);	// Delete the old root.
		}
	}

	/**
	 * This function removes an orphan node. An orphan is defined as a leaf that results from
	 * a network node in the network during the construction of the tree.
	 */
	private void removeOrphan(TMutableNode orphan)
	{
		assert(orphan.isLeaf());

		// Remove this bad node from the tree.
		TMutableNode parent = orphan.getParent();
		parent.removeChild(orphan, false);

		// In case parent has now one child, suppress it.
		if (parent.getChildCount() == 1) {
			suppressBinaryNode(parent);
		}
	}

	/**
	 * This function returns an iterable list of network clusters of this tree.
	 */
	public Iterable<NetworkCluster<T>> getClusters()
	{
		if (_clusters == null) {
			generateClusters();
		}

		return _clusters;
	}


    public Map<NetNode<T>, NetworkCluster<T>> generateNodeClusters()
    {
        Map<NetNode<T>, NetworkCluster<T>> map = new Hashtable<NetNode<T>, NetworkCluster<T>>();

        for (NetNode<T> node : _net.bfs()) {
            getNodeCluster(node, map);
        }

        return map;
    }

	/**
	 * This function generates all clusters for this tree. Trivial clusters, i.e. clusters for single-leaves
	 * are not included.
	 */
	private void generateClusters()
	{
		Map<NetNode<T>, NetworkCluster<T>> map = new Hashtable<NetNode<T>, NetworkCluster<T>>();

		_clusters = new LinkedList<NetworkCluster<T>>();
		for (NetNode<T> node : _net.bfs()) {
			// Compute the cluster for this node.
			NetworkCluster<T> nc = getNodeCluster(node, map);

			// Add it to _clusters.
			if (!nc.isEmpty() && !nc.isTrivial() && !_clusters.contains(nc)) {
				_clusters.add(nc);
			}
		}
	}

	/**
	 * This function returns a cluster under <code>node</code> in this tree.
	 *
	 * @param node: The root of the subtree we want to find its cluster.
	 * @param map: An entry for <code>node</code> is updated after its cluster is found.
	 * The purpose of a map is to avoid re-computation of clusters for the same node.
	 */
	private NetworkCluster<T> getNodeCluster(NetNode<T> node, Map<NetNode<T>, NetworkCluster<T>> map)
	{
		NetworkCluster<T> nc = map.get(node);
		if (nc == null) {
			// There's no entry for this node. Compute a cluster for it.
			if (node.isLeaf()) {
				List<NetNode<T>> leaves = new LinkedList<NetNode<T>>();

				leaves.add(node);
				nc = new NetworkCluster<T>(_net, leaves);	// Create a single-leaf cluster.
			}
			else {
				nc = new NetworkCluster<T>(_net, null);

				// Take the union of clusters of its children.
				for (NetNode<T> child : node.getChildren()) {
					if (child.isTreeNode()) {
						nc.union(getNodeCluster(child, map));
					}
					else {
						NetworkEdge<T> ne = new NetworkEdge<T>(node, child);
						if (_edges.contains(ne)) {
							nc.union(getNodeCluster(child, map));	// Only combine with network node in the tree.
						}
					}
				}
			}

			// Put an entry for this node in the map to avoid re-computation.
			map.put(node, nc);
		}

		return nc;
	}

	/**
	 * This method checks if two trees generated by two different sets of network edges are equal.
	 *
	 * @param nt: The tree to be compared with this tree.
	 *
	 *  @return <code>true</code> if the two trees are equal; <code>false</code> otherwise.
	 */
	public boolean equals(Object nt)
	{
		assert(nt instanceof NetworkTree);

		// Generate clusters for comparing the two trees, if they are not generated yet.
		if (_clusters == null) {
			generateClusters();
		}
		if(((NetworkTree) nt)._clusters == null) {
			((NetworkTree) nt).generateClusters();
		}

		// The two trees are equal if they have the exactly identical set of clusters.
		if (_clusters.size() != ((NetworkTree) nt)._clusters.size()) {
			return false;
		}
		else {
			boolean identical = true;
			Iterator<NetworkCluster<T>> it = _clusters.iterator();

			while (it.hasNext() && identical) {
				NetworkCluster<T> nc = it.next();
				if (!((NetworkTree) nt)._clusters.contains(nc)) {
					identical = false;
				}
			}

			return identical;
		}
	}

	// Data members
	private Network<T> _net;				// The network to build _tree.
	private List<NetworkEdge<T>> _edges;	// List of network edges used to build _tree.
	private List<NetworkCluster<T>> _clusters;	// List of clusters for this tree.
}
