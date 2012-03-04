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

package edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent;

import edu.rice.cs.bioinfo.programs.phylonet.algos.lca.SchieberVishkinLCA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/2/11
 * Time: 3:52 PM
 * To change this template use File | Settings | File Templates.
 */
public class CoalescentCounter
{
    /**
	 * Used internally for storing an edge/cluster pair.
	 *
	 * @author Derek Ruths
	 */
	protected class CEPair {

		public TNode cluster;
		public TNode edge;

		public CEPair(TNode edge, TNode cluster) {
			this.edge = edge;
			this.cluster = cluster;
		}

		public void set(TNode edge, TNode cluster) {
			this.edge = edge;
			this.cluster = cluster;
		}

		public int hashCode() {
			return edge.getID();
		}

		public boolean equals(Object o) {
			if(!(o instanceof CEPair)) {
				return false;
			}

			CEPair p2 = (CEPair) o;

			return (cluster == p2.cluster) && (edge == p2.edge);
		}
	}

	/**
	 * This method counts the number of coalescents that are congruent with the two trees
	 * specified.
	 */
	public int countCoalescents(Tree species_tree, Tree gene_tree) {

		/*
		 * This implementation makes use of the fact that a cluster in a tree is actually the same as a
		 * clade containing two or more leaves.  As a result, nowhere do we explicitly create clusters.
		 * Instead, we take the root of the clade to be the handle to the cluster.
		 */

		SchieberVishkinLCA lca_st = new SchieberVishkinLCA(species_tree);

		Hashtable<CEPair,Integer> ro = new Hashtable<CEPair,Integer>();

		List<TNode> lowest_edges = new LinkedList<TNode>();

		/*
		 * Compute the lowest edges in the tree
		 */
		computeLowestEdges(lca_st, gene_tree.getRoot(), lowest_edges);

		/*
		 * Compute ro for all the edge/cluster pairs.
		 */
		computeRo(lca_st, gene_tree.getRoot(), lowest_edges, ro);

		return ro.get(new CEPair(species_tree.getRoot(), gene_tree.getRoot()));
	}

	/**
	 * Compute the lowest edges in the tree that can contain coalescent events.
	 *
	 * @param lca_st
	 * @param gt_cluster
	 * @param ledges
	 */
	protected void computeLowestEdges(SchieberVishkinLCA lca_st, TNode gt_cluster, List<TNode> ledges) {

		if(gt_cluster.isLeaf()) {
			return;
		}

		TNode lca = getClusterLCA(lca_st, gt_cluster);

		if(ledges.size() == 0) {
			ledges.add(lca);
		} else {
			boolean is_less = true;
			Iterator<TNode> eit = ledges.iterator();
			while(eit.hasNext()) {
				TNode n = eit.next();

				if(!disjoint_paths(lca,n)) {
					if(leq(lca, n)) {
						is_less = true;
						eit.remove();
						break;
					} else {
						is_less = false;
						break;
					}
				}
			}

			if(is_less) {
				ledges.add(lca);
			}
		}

		// recur to the children
		for(TNode child : gt_cluster.getChildren()) {
			computeLowestEdges(lca_st, child, ledges);
		}

	}

	/**
	 * @return the lowest common ancestor in the species tree of the gene tree cluster.
	 */
	protected TNode getClusterLCA(SchieberVishkinLCA lca_st, TNode gt_cluster) {
		List<TNode> leaves = new LinkedList<TNode>();
		getLeaves(gt_cluster, leaves);

		Iterator<TNode> lit = leaves.iterator();

		TNode lca = lca_st.getTree().getNode(lit.next().getName());
		while(lit.hasNext()) {
			lca = lca_st.getLCA(lca, lca_st.getTree().getNode(lit.next().getName()));
		}

		return lca;
	}

	/**
	 * Recursively compute ro for every edge/cluster pair.  This is implemented using head recursion.
	 */
	protected void computeRo(SchieberVishkinLCA lca_st, TNode gt_cluster, List<TNode> lowest_edges, Hashtable<CEPair,Integer> ro) {

		if(gt_cluster.isLeaf()) {
			return;
		}

		Iterator<TNode> pc_it = get_Pc(lca_st, gt_cluster);

		// if this is the lowest cluster
		if(gt_cluster.getLeafCount() <= 2) {

			while(pc_it.hasNext()) {
				TNode edge = pc_it.next();
				ro.put(new CEPair(edge,gt_cluster),1);
			}
		} else {

			// compute children first
			for(TNode child : gt_cluster.getChildren()) {
				computeRo(lca_st, child, lowest_edges, ro);
			}

			// compute this cluster's score
			while(pc_it.hasNext()) {
				TNode edge = pc_it.next();

				if(lowest_edges.contains(edge)) {
					ro.put(new CEPair(edge,gt_cluster), 1);
				} else {
					// compute the product of the sums
					CEPair cp = new CEPair(null,null);
					int sum_prod = 0;

					for(TNode child : gt_cluster.getChildren()) {

						// skip leaves since they aren't clusters
						if(child.isLeaf()) {
							continue;
						}

						int sum = 0;

						Iterator<TNode> child_pc_it = get_Pc(lca_st, child);

						while(child_pc_it.hasNext()) {
							TNode cedge = child_pc_it.next();

							if(leq(cedge,edge)) {
								cp.set(cedge,child);
								sum += ro.get(cp);
							}
						}

						if(sum_prod == 0) {
							sum_prod = sum;
						} else {
							sum_prod *= sum;
						}
					}

					ro.put(new CEPair(edge,gt_cluster), sum_prod);
				}
			}
		}
	}

	/**
	 * @return do nodes <code>e1</code> and <code>e2</code> belong to disjoint paths to the root?
	 */
	protected static boolean disjoint_paths(TNode e1, TNode e2) {

		TNode e1t = e1;

		while(e1 != null) {
			if(e1 == e2) {
				return false;
			}

			e1 = e1.getParent();
		}

		e1 = e1t;

		while(e2 != null) {
			if(e1 == e2) {
				return false;
			}

			e2 = e2.getParent();
		}

		return true;
	}

	/**
	 * @return is node <code>e1</code> less than or equal to node <code>e2</code>?  Here "less than or equal to"
	 * refers to whether e1 is either
	 * <UL>
	 * 	<LI> e1's path to the root contains e2
	 *  <LI> e1's path to the root does not pass through e2 AND e2's path to the root does not pass through e1.
	 * </UL>
	 */
	protected static boolean leq(TNode e1, TNode e2) {

		TNode e1t = e1;

		while(e1 != null) {
			if(e1 == e2) {
				return true;
			}

			e1 = e1.getParent();
		}

		e1 = e1t;

		while(e2 != null) {
			if(e1 == e2) {
				return false;
			}

			e2 = e2.getParent();
		}

		return true;
	}

	/**
	 * Return the edges in the species tree that the gene tree cluster could have coalescent in.
	 */
	protected Iterator<TNode> get_Pc(SchieberVishkinLCA lca_st, TNode gt_cluster) {

		if(gt_cluster.isLeaf()) {
			throw new RuntimeException("get_Pc should not be called on a leaf!");
		}

		final TNode flca = getClusterLCA(lca_st, gt_cluster);

		return new Iterator<TNode>() {

			protected TNode next = flca;

			public boolean hasNext() {
				return (next != null);
			}

			public TNode next() {
				TNode result = next;
				next = next.getParent();

				return result;
			}

			public void remove() {
				throw new RuntimeException("Remove not suppoted");
			}
		};
	}

	/**
	 * Put the leaves under <code>node</code> into the list.
	 */
	protected void getLeaves(TNode node, List<TNode> leaves) {

		if(node.isLeaf()) {
			leaves.add(node);
		} else {
			for(TNode child : node.getChildren()) {
				getLeaves(child, leaves);
			}
		}
	}

}
