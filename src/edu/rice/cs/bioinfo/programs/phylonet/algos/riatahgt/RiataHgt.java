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

package edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt;

import edu.rice.cs.bioinfo.programs.phylonet.algos.lca.SchieberVishkinLCA;
import edu.rice.cs.bioinfo.programs.phylonet.algos.mast.SteelWarnowMAST;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.Serializable;
import java.util.*;

/**
 * This class is a wrapper class for the two classes SingleRiataHgt and MultipleRiataHgt,
 * which computes single and multiple HGT scenarios.
 * The purpose of this class is to provide a consistent interface for computing HGT
 * scenarios between two trees.
 */
public class RiataHgt implements Serializable {


	/**
	 * This function constructs an object of RiataHgt for computing HGT events.
	 * Flags _collapsed and _refined are set to <code>true</code> by default. If the user does not want
	 * to collapse and/or refine the trees before computing HGT events, call the function disableCollapse()
	 * and/or disableRefine(), respectively.
	 */
	public RiataHgt()
	{
		_collapsed = true;
		_refined = true;
	}

	/**
	 * This function tells RiataHgt to collapse the trees before computing HGT events.
	 */
	public void enableCollapse()
	{
		_collapsed = true;
	}

	/**
	 * This function tells RiataHgt not to collapse the trees before computing HGT events.
	 */
	public void disableCollapse()
	{
		_collapsed = false;
	}

	/**
	 * This function tells RiataHgt to refine the trees before computing HGT events.
	 */
	public void enableRefine()
	{
		_refined = true;
	}

	/**
	 * This function tells RiataHgt not to refine the trees before computing HGT events.
	 */
	public void disableRefine()
	{
		_refined = false;
	}

	/**
	 * This function sets a string used to label internal node. Its default value is "I".
	 *
	 * @param pfx: The string used as a prefix for internal nodes.
	 */
	public void setPrefix(String pfx)
	{
		if (pfx != "") {
			_prefix = pfx;
		}
	}

	/**
	 * This function computes all HGT scenarios for the two species and gene trees. This function should
	 * be called before the functions getMinimumSolutionSize and/or getMinHgtScenarios can be called.
	 *
	 * @param species, gene: The species and gene trees.
	 */
	public void computeHgt(Tree species, Tree gene)
	{
		// Preprocess the species and gene trees before computing HGT events.
		preprocessTrees(species, gene);

		// Compute all possible HGT events.
		_result_tree = new STITree<List<HgtScenario>>();
		if (_result_tree.isEmpty()) {
			_result_tree.createRoot();
		}

		_result_tree.getRoot().setName(species.getRoot().getName());
		computeMultipleHgtHelper(_species_copy, _gene_copy, _result_tree.getRoot());
	}

	public void computeSingleHgt(Tree species, Tree gene)
	{
		// Preprocess the species and gene trees before computing HGT events.
		preprocessTrees(species, gene);

		// Compute all possible HGT events.
		_result_tree = new STITree<List<HgtScenario>>();
		if (_result_tree.isEmpty()) {
			_result_tree.createRoot();
		}

		_result_tree.getRoot().setName(species.getRoot().getName());
		computeSingleHgtHelper(_species_copy, _gene_copy, _result_tree.getRoot());
	}

	private void computeMultipleHgtHelper(MutableTree st, MutableTree gt, STINode<List<HgtScenario>> node)
	{
		MultipleRiataHgt mrh = new MultipleRiataHgt();
		List<MutableTree> subst_copies = new LinkedList<MutableTree>();
		List<MutableTree> subgt_copies = new LinkedList<MutableTree>();

		// Collapse st and gt first.
		if (_collapsed) {
			collapseClusters(st, gt, subst_copies, subgt_copies);
		}

		// Compute HGT events for st and gt.
		mrh.computeMultipleHgt(st, gt);
		node.setData(mrh.getMinHgtScenarios());

		// Compute for collapsed subtrees in st and gt.
		if (_collapsed) {
			for (int i = 0; i < subst_copies.size(); i++) {
				String name = subst_copies.get(i).getRoot().getName();
				STINode<List<HgtScenario>> child = node.createChild(name);

				computeMultipleHgtHelper(subst_copies.get(i), subgt_copies.get(i), child);
			}
		}
	}

	private void computeSingleHgtHelper(MutableTree st, MutableTree gt, STINode<List<HgtScenario>> node)
	{
		SingleRiataHgt srh = new SingleRiataHgt();
		List<MutableTree> subst_copies = new LinkedList<MutableTree>();
		List<MutableTree> subgt_copies = new LinkedList<MutableTree>();

		// Collapse st and gt first.
		collapseClusters(st, gt, subst_copies, subgt_copies);

		// Compute HGT events for st and gt.
		srh.computeHGT(st, gt);
		List<HgtScenario> temp = new LinkedList<HgtScenario>();
		HgtScenario scenario = new HgtScenario();
		for (HgtEvent event : srh.getHGTEvents()) {
			scenario.addEvent(event);
		}
		temp.add(scenario);

		node.setData(temp);

		// Compute for collapsed subtrees in st and gt.
		for (int i = 0; i < subst_copies.size(); i++) {
			String name = subst_copies.get(i).getRoot().getName();
			STINode<List<HgtScenario>> child = node.createChild(name);

			computeSingleHgtHelper(subst_copies.get(i), subgt_copies.get(i), child);
		}
	}

	/**
	 * Enumerate all solutions computed by RiataHgt.
	 *
	 * @return: A list of HgtScenarios computed by RiataHgt.
	 */
	public List<HgtScenario> enumerateSolutions() {
		return enumerateSubsolutions(_result_tree.getRoot());
	}

	/**
	 * Get a (strict) consensus network from HGT events of different gene trees.
	 *
	 * @param allGenesEvents: Each member of the list is a set HgtEvents that appear
	 * @return: A scenario that consists of events that appear in all scenarios.
	 */
	public static List<HgtEvent> getConsensusNetwork(List<List<HgtEvent>> allGenesEvents) {
		List<HgtEvent> result = new LinkedList<HgtEvent>();

		if (allGenesEvents.size() == 0) {
			return result;
		}
		else if (allGenesEvents.size() == 1) {
			return allGenesEvents.get(0);
		}

		for (HgtEvent event : allGenesEvents.get(0)) {
			boolean all = true;
			for (int i = 1; i < allGenesEvents.size(); i++) {
				if (!allGenesEvents.get(i).contains(event)) {
					all = false;
					break;
				}
			}

			if (all) {
				result.add(event);
			}
		}

		return result;
	}

	/**
	 * A helper function for enumerateSubsolutions. This function is implemented as a recursive procedure.
	 *
	 * @param node: The node we want to get the subsolutions.
	 * @return: The list subsolutions under this node.
	 */
	private List<HgtScenario> enumerateSubsolutions(STINode<List<HgtScenario>> node) {
		if (node.isLeaf()) {
			return node.getData();
		}
		else {
			List<HgtScenario> scenarios = new LinkedList<HgtScenario>();

			for (HgtScenario hs : node.getData()) {
				scenarios.add(hs);
			}

			for (STINode<List<HgtScenario>> child : node.getChildren()) {
				if (scenarios.size() == 0) {
					scenarios = enumerateSubsolutions(child);
				}
				else {
					List<HgtScenario> temp = new LinkedList<HgtScenario>();
					List<HgtScenario> childScenarios = enumerateSubsolutions(child);
					for (HgtScenario hs1 : scenarios) {
						for (HgtScenario hs2 : childScenarios) {
							HgtScenario hs = new HgtScenario();

							hs.addEvents(hs1);
							hs.addEvents(hs2);
							temp.add(hs);
						}
					}

					scenarios = temp;
				}
			}

			return scenarios;
		}
	}

	/**
	 * This function returns the size of a minium HGT scenario (i.e., the minimum number of HGT events that
	 * Riata can detect to induce the gene tree from the species tree). This function can only be called after
	 * the function computeHgt is called.
	 *
	 * @return The size of a minimum HGT scenario.
	 */
	public int getMinimumSolutionSize()
	{
		int count = 0;

		for (STINode<List<HgtScenario>> node : _result_tree.getNodes()) {
			if (!node.getData().isEmpty()) {
				HgtScenario hs = node.getData().get(0);
				count += hs.getEvents().size() - hs.countBadEvents();
			}
		}

		return count;
	}

	/**
	 * This function returns the number of minimal solutions.
	 */
	public int getNumberOfMinSolutions()
	{
		STINode<List<HgtScenario>> root = _result_tree.getRoot();

		if (root.getData().isEmpty() && root.getChildCount() == 0) {
			return 0;	// Two trees are identical.
		}
		else {
			int count = 1;

			for (STINode<List<HgtScenario>> node : _result_tree.getNodes()) {
				count *= getNodeSolutionCount(node);
			}

			return count;
		}
	}

	/**
	 * This function counts the number of (partial) solutions within each node.
	 */
	private int getNodeSolutionCount(STINode<List<HgtScenario>> node)
	{
		if (node.getData().isEmpty()) {
			return 1;
		}
		else {
			int count = 0;

			for (HgtScenario hs : node.getData()) {
				int temp = 1;

				for (HgtEvent event : hs.getEvents()) {
					TNode pnode = event._primary_node;

					if (pnode != null && !_counted_primary_nodes.contains(pnode)) {
						temp *= event.getNumberOfAlternatives();
						_counted_primary_nodes.add(pnode);
					}
				}

				count += temp;
			}

			return count;
		}
	}

	/**
	 * Computes the weights of events returned by computeHgt. The weight of an event
	 * is defined as the number of solutions in which this event appear. This weight
	 * might need to be normalized by the total number of solutions returned.
	 */
	public Map<HgtEvent, Integer> computeEventWeights() {
		Map<HgtEvent, Integer> weights = new Hashtable<HgtEvent, Integer>();
		int solutionCount = getNumberOfMinSolutions();

		for (STINode<List<HgtScenario>> node : _result_tree.getNodes()) {
			if (!node.getData().isEmpty()) {
				int size = node.getData().size();	// The number of alternate subsolutions for this node.
				for (HgtScenario hs : node.getData()) {
					for (HgtEvent event : hs.getEvents()) {
						int w = solutionCount / size;	// The number of solutions this event will appear in.
						if (!weights.containsKey(event)) {
							weights.put(event, w);
						}
						else {
							w += weights.get(event);	// Increase the weight for this event.
							weights.remove(event);
							weights.put(event, w);
						}
					}
				}
			}
		}

		return weights;
	}

	/**
	 * This function returns the solution tree that contains HGT scenarios for a pair of species and gene trees.
	 * Each node in the tree corresponds to an independent component. Therefore, to obtain a complete set of
	 * HGT events for the whole species and gene trees, select a scenario from each node.
	 */
	public STITree<List<HgtScenario>> getSolutionTree()
	{
		return _result_tree;
	}

	/**
	 * Returns the original, unmodified (except for labeling) species tree used
	 * by the last RIATA-HGT calculation.
	 *
	 * @return	the last species tree used
	 */
	public Tree getSpeciesTree() {
		return _species_original;
	}

	/**
	 * This function prints the solution.
	 */
	public String toString()
	{
		String output = "There are " + getComponentCount() + " component(s), which account(s) for " + getNumberOfMinSolutions() + " solution(s), each of size " + getMinimumSolutionSize() + "\n";
		output += toString(_result_tree.getRoot(), 0);

		return output;
	}

	public int getComponentCount() {
		int component_count = 0;
		for (STINode<List<HgtScenario>> node : _result_tree.getNodes()) {
			if (!node.getData().isEmpty()) {
				component_count++;
			}
		}
		return component_count;
	}

	private String toString(STINode<List<HgtScenario>> node, int tab_level)
	{
		String output = "";

		if (!node.getData().isEmpty()) {
			output = "-----------------------------------------------------------------------------------------------------";
			output += "\n" + printTabs(tab_level) + "Component " + node.getName() + ":\n";
			int i = 1;

			for (HgtScenario hs : node.getData()) {
				output += printTabs(tab_level) + "Subsolution" + i + ":\n";
				output += hs.toString(tab_level + 1);

				i++;
			}
		}

		for (STINode<List<HgtScenario>> child : node.getChildren()) {
			output += toString(child, tab_level + 1);
		}

		return output;
	}

	private String printTabs(int tab_level) {
		String tabs = "";
		for (int i = 0; i < tab_level; i++) {
			tabs += "\t";
		}

		return tabs;
	}

	/**
	 * This function preprocesses the two trees before RiataHgt computes HGT events.
	 *
	 * @param species, gene: The two trees to be preprocessed.
	 */
	public void preprocessTrees(Tree species, Tree gene)
	{
		// Refine the two trees, if the flag _refined is enabled.
		ExMultipleMasts emm = new ExMultipleMasts();

		if (_refined) {
			emm.refine(species, gene);
		}

		_species_original = species;
		Trees.autoLabelNodes((MutableTree) _species_original, _prefix);
		_species_copy = new STITree(species);

		// Clear names of all internal nodes, preapring for collapsing two trees.
		for (TNode node : gene.getNodes()) {
			if (!node.isLeaf()) {
				((TMutableNode) node).setName(TNode.NO_NAME);
			}
		}
		_gene_copy = new STITree(gene);

		if (_collapsed) {
			emm.collapse(_species_copy, _gene_copy);
		}

		Trees.removeBinaryNodes(_species_copy);
		Trees.removeBinaryNodes(_gene_copy);

		// Collapse chains.
		List<Chain> stChains = new LinkedList<Chain>();
		List<Chain> gtChains = new LinkedList<Chain>();
		findChains((MutableTree) _species_copy, (MutableTree) _gene_copy, stChains, gtChains);
		collapseChains((MutableTree) _species_copy, (MutableTree) _gene_copy, stChains, gtChains);

		// Relabels nodes in the gene according to nodes in _gene_copy.
		for (TNode it : _gene_copy.getNodes()) {
			TNode peer = gene.getNode(it.getID());
			if (peer != null && peer.getName().equals(peer.NO_NAME)) {
				((STINode) peer).setName(it.getName());
			}
		}
	}

	/**
	 * This function is for collapsing identical clades between _species_copy and _gene_copy. This function
	 * assumes that the function preprocessTrees() has been called before it is called, so that all internal nodes
	 * in _species_tree have a name.
	 */
	private void collapseClusters(MutableTree st, MutableTree gt, List<MutableTree> subst_copies, List<MutableTree> subgt_copies)
	{
		TNode old = null;
		LinkedList<TNode> temp = new LinkedList<TNode>();
		for (TNode child : st.getRoot().getChildren()) {
			temp.add(child);
		}

		int size = temp.size();
		for (int i = 0; i < size; i++) {
			old = temp.get(i);
			collapseNode(old, gt, subst_copies, subgt_copies);
		}
	}

	/**
	 * This function collapse the subtree rooted at <code>node</code> of the species tree. If the clade under
	 * <code>node</code> is not identical to any clade in the gene tree, then the function recursively collapes
	 * its children.
	 */
	private void collapseNode(TNode node, MutableTree gt, List<MutableTree> subst_copies, List<MutableTree> subgt_copies)
	{
		if (node.isLeaf()) {
			return;
		}

		MutableTree sst = new STITree(node);	// The subtree rooted at node on the species tree.

		// Find the LCA in the gene tree for the set of leaves in sst.
		Set<TNode> gst_leaves = new HashSet<TNode>();
		for (TNode n : sst.getNodes()) {
			if (n.isLeaf()) {
				gst_leaves.add(gt.getNode(n.getName()));
			}
		}

		SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(gt);
		TNode lca = lcaFinder.getLCA(gst_leaves);

		// Check if we can collapse node and lca?
		MutableTree gst = new STITree(lca);
		if (sst.getLeafCount() == gst.getLeafCount()) {	// Two clades are equal; collapse.
			subst_copies.add(sst);
			subgt_copies.add(gst);

			// Remove node.
			String name = node.getName();
			TMutableNode parent = (TMutableNode) node.getParent();

			parent.removeChild((TMutableNode) node, false);
			parent.createChild(name);

			// Remove lca.
			parent = (TMutableNode) lca.getParent();
			parent.removeChild((TMutableNode) lca, false);
			parent.createChild(name);
		}
		else {	// Two clades are not equal; collapse their children.
			TNode old = null;

			LinkedList<TNode> temp = new LinkedList<TNode>();

			for (TNode child : node.getChildren()) {
				temp.add(child);
			}

			for (int i = 0; i < temp.size(); i++) {
				old = temp.get(i);
				collapseNode(old, gt, subst_copies, subgt_copies);
			}
		}
	}

	/**
	 * Finds maximal identical chains in st and gt.
	 *
	 * @param: st, gt: The two trees where we want to find chains.
	 * @param: stChains, gtChains: Two lists of chains found in st and gt.
	 *
	 * ASSUMPTION: Collapsing identical clades must be done before this function
	 * is called.
	 */
	public void findChains(MutableTree st, MutableTree gt, List<Chain> stChains, List<Chain> gtChains)
	{
		List<TNode> startNodes = new LinkedList<TNode>();
		List<Integer> startNodeDepths = new LinkedList<Integer>();
		Set<TNode> checkedNodes = new HashSet<TNode>();

		// Get leaves as the starting set.
		for (TNode node : st.getNodes()) {
			if (node.isLeaf()) {
				int depth = getNodeDepth(node);
				addToStartNodes(startNodes, startNodeDepths, node, depth);
			}
		}

		// Begin to find maximal chains.
		while (!startNodes.isEmpty()) {
			TNode v = startNodes.remove(0);
			int v_depth = startNodeDepths.remove(0);

			if (checkedNodes.contains(v) || v.isRoot()) {
				continue;
			}
			else {
				checkedNodes.add(v);
			}

			TNode u = v.getParent();
			int u_depth = v_depth - 1;
			Set<String> L_uv = new HashSet<String>();	// Holds L(u) - L(v).

			if (checkedNodes.contains(u)) {
				continue;	// different maximal chains do not share nodes.
			}

			// First, only consider the cases where all children of u, except v, are leaves.
			boolean on_chain = true;
			for (TNode child : u.getChildren()) {
				if (child != v) {
					if (child.isLeaf()) {
						L_uv.add(child.getName());
					}
					else {	// u cannot be on a chain.
						// but u can be a cut on a chain.
						addToStartNodes(startNodes, startNodeDepths, u, u_depth);
						on_chain = false;
						break;
					}
				}
			}

			if (on_chain) {
				// L_uv now equals L(u) - L(v), find edge e' in gt.
				TNode v_prime = findStartEdge(gt, L_uv);
				if (v_prime == null) {
					// No such edge exists, continue with another node.
					addToStartNodes(startNodes, startNodeDepths, u, u_depth);
					continue;
				}
				else {
					// There's such edge e', but need to check topology.
					TNode u_prime = v_prime.getParent();
					boolean agreed = true;
					for (TNode child : u_prime.getChildren()) {
						if (child != v_prime) {
							if (!child.isLeaf()) {
								agreed = false;
								break;
							}
						}
					}

					if (!agreed) { // e and e' cannot be extended to become chains.
						addToStartNodes(startNodes, startNodeDepths, u, u_depth);
						continue;
					}
					else {
						// Try to extend e and e'.
						Chain st_chain = new Chain();
						Chain gt_chain = new Chain();

						st_chain._cut = v;
						st_chain._bottom = st_chain._top = u;
						gt_chain._cut = v_prime;
						gt_chain._bottom = gt_chain._top = u_prime;

						extendChains(st_chain, gt_chain, checkedNodes);

						// Remember the chains for later collapsing.
						if (st_chain.getLength() >= 1) {
							stChains.add(st_chain);
							gtChains.add(gt_chain);

							checkedNodes.add(u);
						}
						else {
							// u can be v in another chain.
							addToStartNodes(startNodes, startNodeDepths, u, u_depth);
						}
					}
				}
			}
		}
	}

	/**
	 * Replace identical chains found by findChains by abc-chains.
	 *
	 * @param st, gt: The species and gene trees.
	 * @param stChains, gtChains: Chains in the species and gene trees.
	 */
	public void collapseChains(MutableTree st, MutableTree gt, List<Chain> stChains, List<Chain> gtChains)
	{
		if (stChains.isEmpty()) {
			return;
		}

		for (int i = 0; i < stChains.size(); i++) {
			Chain st_ch = stChains.get(i);
			Chain gt_ch = gtChains.get(i);

			if (st_ch.getLength() >= 3) {
				// Collapse in the species tree.
				TNode v = st_ch._cut;
				TMutableNode u = (TMutableNode) st_ch._bottom;

				for (int j = 0; j < 3; j++) {
					List<TNode> children = new LinkedList<TNode>();
					for (TNode child : u.getChildren()) {
						if (child != v) {
							children.add(child);
						}
					}

					if (children.size() > 1) {
						((TMutableNode) children.get(0)).setName(_chain_prefix + i + "_" + j);
						for (int k = 1; k < children.size(); k++) {
							u.removeChild((TMutableNode) children.get(k), false);
						}
					}

					v = u;
					u = u.getParent();
				}

				TMutableNode topChain = (TMutableNode) v;

				while (v != st_ch._top) {
					List<TNode> children = new LinkedList<TNode>();
					for (TNode child : u.getChildren()) {
						if (child != v) {
							children.add(child);
						}
					}

					for (int k = 0; k < children.size(); k++) {
						u.removeChild((TMutableNode) children.get(k), false);
					}

					v = u;
					u = u.getParent();
				}

				// Set v's name to the top abc-chain node. VERY IMPORTANT.
				String topName = st_ch._top.getName();
				((TMutableNode) st_ch._top).setName(STITree.NO_NAME);
				topChain.setName(topName);

				// Collapse in the gene tree.
				v = gt_ch._cut;
				u = (TMutableNode) gt_ch._bottom;

				for (int j = 0; j < 3; j++) {
					List<TNode> children = new LinkedList<TNode>();
					for (TNode child : u.getChildren()) {
						if (child != v) {
							children.add(child);
						}
					}

					if (children.size() > 1) {
						((TMutableNode) children.get(0)).setName(_chain_prefix + i + "_" + j);
						for (int k = 1; k < children.size(); k++) {
							u.removeChild((TMutableNode) children.get(k), false);
						}
					}

					v = u;
					u = u.getParent();
				}

				topChain = (TMutableNode) v;

				while (v != gt_ch._top) {
					List<TNode> children = new LinkedList<TNode>();
					for (TNode child : u.getChildren()) {
						if (child != v) {
							children.add(child);
						}
					}

					for (int k = 0; k < children.size(); k++) {
						u.removeChild((TMutableNode) children.get(k), false);
					}

					v = u;
					u = u.getParent();
				}

				// Set v's name to the top abc-chain node. VERY IMPORTANT.
				topName = gt_ch._top.getName();
				((TMutableNode) gt_ch._top).setName(STITree.NO_NAME);
				topChain.setName(topName);
			}
		}

		// Remove any binary nodes as results of collapsing chains.
		Trees.removeBinaryNodes(st);
		Trees.removeBinaryNodes(gt);
	}

	/**
	 * Returns the set of leaves of the subtree rooted at subroot.
	 *
	 * @param subroot: root for the subtree.
	 * @return: set of leaves under subroot.
	 */
	private Set<String> getSubtreeLeaves(TNode subroot)
	{
		Set<String> ret = new HashSet<String>();

		if (subroot.isLeaf()) {
			ret.add(subroot.getName());
			return ret;
		}
		else {
			for (TNode child : subroot.getChildren()) {
				ret.addAll(getSubtreeLeaves(child));
			}
			return ret;
		}
	}

	/**
	 * Find the corresponding edges e' = (u', v') such that L(u') - L(v') = leaves.
	 *
	 * @param gt: The tree
	 * @return v' of e' (so that e' can be found by taking u' = parent(v').
	 */
	private TNode findStartEdge(MutableTree gt, Set<String> names)
	{
		SchieberVishkinLCA sv = new SchieberVishkinLCA(gt);

		Set<TNode> leaves = new HashSet<TNode>();
		for (String str : names) {
			TNode n = gt.getNode(str);
			leaves.add(n);
		}

		TNode u = sv.getLCA(leaves);	// u' in GT.
		Set<String> L_u = getSubtreeLeaves(u);

		if (L_u.size() == leaves.size()) {
			u = u.getParent();
			L_u = getSubtreeLeaves(u);
		}

		for (TNode v : u.getChildren()) {
			Set<String> L_v = getSubtreeLeaves(v);
			Set<String> temp = new HashSet<String>(L_u);
			temp.removeAll(L_v);

			// Check names = temp.
			if (temp.containsAll(names) && names.containsAll(temp)) {
				return v;
			}
		}

		return null;
	}

	/**
	 * Extends the chain upward until it cannot be extended.
	 *
	 * @param st_chain, gt_chain: Two initial chains to be extended.
	 * @param checkedNodes: To update nodes that we don't need to look at as a starting node
	 * of a chain.
	 */
	private void extendChains(Chain st_chain, Chain gt_chain, Set<TNode> checkedNodes)
	{
		while (true) {
			// Reach root.
			if (st_chain._top.isRoot() || gt_chain._top.isRoot()) {
				return;
			}

			// Else, extend.
			Set<String> u_siblings = new HashSet<String>();
			for (TNode child : st_chain._top.getParent().getChildren()) {
				if (child != st_chain._top) {
					if (child.isLeaf()) {
						u_siblings.add(child.getName());
					}
					else {
						// Cannot extend.
						return;
					}
				}
			}

			Set<String> u_prime_siblings = new HashSet<String>();
			for (TNode child : gt_chain._top.getParent().getChildren()) {
				if (child != gt_chain._top) {
					if (child.isLeaf()) {
						u_prime_siblings.add(child.getName());
					}
					else {
						// Cannot extend.
						return;
					}
				}
			}

			// Now, go up one node if these siblings agree.
			if (u_siblings.containsAll(u_prime_siblings) && u_prime_siblings.containsAll(u_siblings)) {
				st_chain._top = st_chain._top.getParent();
				gt_chain._top = gt_chain._top.getParent();

				checkedNodes.add(st_chain._top);
			}
			else {
				return;
			}
		}
	}

	/**
	 * Returns the depth of a node.
	 *
	 * @param node: Node to compute its depth.
	 */
	private int getNodeDepth(TNode node)
	{
		int depth = 0;
		while (!node.isRoot()) {
			node = node.getParent();
			depth++;
		}

		return depth;
	}

	/**
	 * Add a new node to startNodes so that the list still maintains
	 * nodes in decreasing depths.
	 *
	 * @param startNodes: List of nodes in decreasing depths
	 * @param depthNodes: Depths of nodes corresponding to startNodes.
	 * @param n: node to be added.
	 * @param d: n's depth.
	 */
	private void addToStartNodes(List<TNode> startNodes, List<Integer> nodeDepths, TNode n, int d)
	{
		if (startNodes.isEmpty()) {
			startNodes.add(n);
			nodeDepths.add(d);
		}
		else {
			int i;
			for (i = 0; i < nodeDepths.size(); i++) {
				if (startNodes.get(i) == n) {
					return;
				}

				if (nodeDepths.get(i) < d) {
					break;
				}
			}

			startNodes.add(i, n);
			nodeDepths.add(i, d);
		}
	}

	// Data members.
	STITree<List<HgtScenario>> _result_tree;	// Holds the final scenario(s).

	private boolean _collapsed;	// Indicates if the two trees should be collapsed.
	private boolean _refined;	// Indicates if the two trees should be refined.

	private SingleRiataHgt _rh;				// This data member is for computing single HGT scenario.
	private MultipleRiataHgt _ermh;	// This data member is for computing multiple HGT scenarios.

	private Tree _species_original;		// Stores the original, "unmodified" version of the species tree.
	private MutableTree _species_copy;	// Stores a copy of the species tree for computing HGT events.
	private MutableTree _gene_copy;		// Stores a copy of the gene tree for computing HGT events.
	private List<TNode> _counted_primary_nodes = new LinkedList<TNode>();

	private String _prefix = "I";		// The prefix used to label internal nodes.
	private String _chain_prefix = "C_I";	// Prefix for abc-chains.
}

/**
 * Class represents an identical chain in two trees.
 */
class Chain {
	// Data members
	public TNode _top, _bottom;	// Top and bottom nodes of chain.
	public TNode _cut;			// Cut node at bottom.

	// Return the length of the chain.
	public int getLength()
	{
		TNode n = _bottom;
		int len = 0;

		while (n != _top) {
			n = n.getParent();
			len++;
		}

		return len;
	}
}

/**
 * This class implements the RIATA-HGT algorithm outlined in the paper
 * <B>RIATA-HGT: A Fast and Accurate Heuristic for Reconstructing Horizontal Gene Transfer</B>
 * by Luay Nakhleh, Derek Ruths, and Li-San Wang.
 *
 * @author Derek Ruths
 */
class SingleRiataHgt {
	// constants
	public static final String OUTGROUP_NAME = "__X";

	/**
	 * This function add an outgroup to tree <code>tr</code>.
	 *
	 * @param: <code>tr</code> is a tree we want to add an outgroup to it.
	 */
	public static Tree addOutgroup(Tree tr)
	{
		assert(tr.isRooted());	// The tree must be rooted.

		STITree<Object> xtr = new STITree<Object>();
		xtr.getRoot().createChild(OUTGROUP_NAME);
		xtr.getRoot().createChild(tr.getRoot());

		return xtr;
	}

	/**
	 * This function is to remove the outgroup we alread added to the tree xtr. Its use is in getting
	 * the real MAST between two original trees.
	 *
	 * @param xtr: The tree with an outgroup.
	 * @return: The tree after an outgroup is removed.
	 */
	public static Tree removeOutgroup(Tree xtr) {
		Tree tr = null;
		for (TNode child : xtr.getRoot().getChildren()) {
			if (child.isLeaf() && child.getName().equals(OUTGROUP_NAME)) {
				continue;	// Ignore the outgroup.
			}
			else {
				tr = new STITree<Object>(child);
				break;		// This node corresponds to the actual MAST.
			}
		}

		return tr;
	}

	// This data structure holds a list of decomposition nodes that have other alternatives.
	// Such nodes are in the case 2-b in the function Decompose in the paper RIATA-HGT: A Fast...
	protected List<TNode> _primary_nodes = new LinkedList<TNode>();
	protected List<TNode> _decomp_nodes = new LinkedList<TNode>();

	/**
	 * This function finds the decomposition node for <code>node</code> if any. It returns null
	 * if no such nodes exist for <code>node</code>.
	 *
	 * @param node:
	 * @return the decomposition node for <code>node</code>.
	 */
	protected TNode getDecompNode(TNode node)
	{
		TNode temp = node;
		 while (true) {
			 TNode parent = temp.getParent();
			 if (parent == null || _decomp_nodes.contains(parent)) {
				 return parent;
			 }
			 else {
				 temp = parent;	// Continue to find.
			 }
		 }
	}

	// fields
	protected Hashtable<STITree<Object>,Tree> _backbone_map = new Hashtable<STITree<Object>,Tree>();
	protected Set<HgtEvent> _hgt_events = null;
	protected Tree _org_gene_tree, _org_species_tree;

	// methods
	public void computeHGT(Tree species_tree, Tree gene_tree) {

		if(!species_tree.isRooted() || !gene_tree.isRooted()) {
			throw new RuntimeException("Trees must be rooted");
		}

		_hgt_events = new HashSet<HgtEvent>();
		_org_gene_tree = gene_tree;
		_org_species_tree = species_tree;
		inner_computeHGT(species_tree, gene_tree);
	}

	private void inner_computeHGT(Tree species_tree, Tree gene_tree) {
		// compute the rooted mast (T')
		Tree xst = SingleRiataHgt.addOutgroup(species_tree);
		Tree xgt = SingleRiataHgt.addOutgroup(gene_tree);

		SteelWarnowMAST mast_alg = new SteelWarnowMAST();
		Tree xmt = mast_alg.computeRMAST(xst, xgt);
		Tree mast_tree = SingleRiataHgt.removeOutgroup(xmt);	//Tree mast_tree = mast_alg.computeRMAST(species_tree, gene_tree);

		STITree<Object> sti_st = new STITree<Object>(species_tree.getRoot(), true);

		/**
		 * If the trees don't overlap at all, there's no HGT information encoded here.
		 */
		if(mast_tree == null || mast_tree.getLeafCount() == 0) {
			return;
		}

		// if the mast is the species tree, there's no difference between the trees
		if(mast_tree.getLeafCount() == species_tree.getLeafCount()) {
			return;
		}

		// U1 = ST - T'
		STITree<Object>[] U1_array = (STITree<Object>[]) createNonMASTSubtrees(species_tree, mast_tree);
		LinkedList<STITree<Object>> U1 = new LinkedList<STITree<Object>>();
		for (int i = 0; i < U1_array.length; i++) {
			U1.add(U1_array[i]);
		}

		LinkedList<Set<String>> U1_leaf_names = new LinkedList<Set<String>>();
		for(int i = 0; i < U1.size(); i++) {
			U1_leaf_names.add(getLeafNames(U1.get(i)));
		}

		// U2 = GT - T'
		STITree<Object>[] U2 = (STITree<Object>[]) createNonMASTSubtrees(gene_tree, mast_tree);

		// Decompose the clades in u2
		Hashtable<STITree<Object>,STITree<Object>> V = new Hashtable<STITree<Object>,STITree<Object>>();
		for(int i = 0; i < U2.length; i++) {
			decompose(U1, U1_leaf_names, U2[i], mast_tree, V);
		}

		U2 = (STITree<Object>[]) new STITree[V.size()];
		Iterator<STITree<Object>> iter = V.keySet().iterator();
		for(int i = 0; iter.hasNext(); i++) {
			U2[i] = iter.next();
		}

		// while ...
		while(!V.isEmpty()) {

			// get one element of V and its corresponding element in U1
			Map.Entry<STITree<Object>,STITree<Object>> entry = V.entrySet().iterator().next();
			//STITree u2 = (STITree) entry.getKey();
			STITree<Object> u1 = entry.getValue();

			// Y = {y \in U2 : L(y) \cap L(u_1) \not= \emptyset}
			// Z = {y|(L(y) - L(u_1)) : y \in Y }
			HashSet<STITree<Object>> Y = new HashSet<STITree<Object>>();
			HashSet<STITree<Object>> Z = new HashSet<STITree<Object>>();
			Set<String> u1_leafset = getLeafNames(u1);
			for(int i = 0; i < U2.length; i++) {
				Set<String> u2_leafset = getLeafNames(U2[i]);
				HashSet<String> s = new HashSet<String>(u2_leafset);

				// check the intersection
				s.addAll(u1_leafset);
				if(s.size() < (u1_leafset.size() + u2_leafset.size())) {
					Y.add(U2[i]);
					STITree<Object> ztree = new STITree<Object>(U2[i].getRoot(), true);

					// s = L(y) - L(u_1)
					s.clear();
					s.addAll(u2_leafset);
					s.removeAll(u1_leafset);
					ztree.constrainByLeaves(s);

					if(ztree.getNodeCount() > 0) {
						Z.add(ztree);
					}
				}
			}

			// V = V - Y
			V.keySet().removeAll(Y);

			// X = {u_1|L(y) : y \in Y}
			Set<STITree<Object>> X = new HashSet<STITree<Object>>();
			for(STITree<Object> t : Y) {
				STITree<Object> u1_cpy = new STITree<Object>(u1.getRoot(), true);
				u1_cpy.constrainByLeaves(getLeafNames(t));

				X.add(u1_cpy);
			}

			// foreach y \in Y ...
			for(STITree<Object> y : Y) {
				// find the x \in X s.t. L(x) \cap L(y) \not= \emptyset
				Set<String> y_leaves = getLeafNames(y);
				STITree<Object> x = null;
				Set<String> s = new HashSet<String>();

				for(STITree<Object> tx : X) {
					x = tx;
					Set<String> x_leaves = getLeafNames(x);

					s.clear();
					s.addAll(y_leaves);
					s.removeAll(x_leaves);

					// if L(x) \cap L(y) \not= \emptyset
					if(s.size() < y_leaves.size()) {
						break;
					}
				}

				inner_computeHGT(x,y);

				addSingleHGT(sti_st, gene_tree, y, U2, mast_tree);
			}
		}
	}

	public Set<HgtEvent> getHGTEvents() {
		return new HashSet<HgtEvent>(_hgt_events);
	}

	/**
	 * @return the set of all the names of the leaves in <code>t</code>.
	 */
	protected Set<String> getLeafNames(Tree t) {
		HashSet<String> s = new HashSet<String>();

		for(TNode n : t.getNodes()) {
			if(n.isLeaf()) {
				s.add(n.getName());
			}
		}

		return s;
	}

	/**
	 * For each leaf in the MAST tree, remove its path to the root in the other
	 * tree.  Retain all subtrees formed from the leaf removal process.
	 *
	 * @param tree
	 * @param mast_tree
	 *
	 * @return
	 */
	public STITree<Object>[] createNonMASTSubtrees(Tree tree, Tree mast_tree) {

		Vector<STITree<Object>> trees = new Vector<STITree<Object>>();
		trees.add(new STITree<Object>(tree.getRoot(), true));

		Set<String> leaf_set = getLeafNames(mast_tree);
		for (String leaf_name : leaf_set) {

			// find the tree that the leaf exists in
			STITree<Object> t = null;
			Iterator<STITree<Object>> j = trees.iterator();
			while (j.hasNext()) {
				t = j.next();

				if (t.getNode(leaf_name) != null) {
					j.remove();
					break;
				}
			}

			// if the tree consists of that leaf, just drop the tree
			// and don't try to get other trees in it.
			if (t.getNodeCount() == 1) {
				continue;
			}

			// remove the leaf from the tree, storing all new trees created
			STINode<Object> n = t.getNode(leaf_name);
			while (n != t.getRoot()) {
				STINode<Object> parent = n.getParent();

				for(STINode<Object> child : parent.getChildren()) {
					// kick off the orphaned subtrees
					if(child != n) {
						trees.add(new STITree<Object>(child, true));
					}
				}

				n = parent;
			}
		}

		// remove extraneous nodes from the trees we have made
		STITree[] result = new STITree[trees.size()];
		int idx = 0;
		for(STITree<Object> nt : trees) {
			result[idx] = nt;

			// clean the tree
			Trees.removeBinaryNodes(result[idx]);

			idx++;
		}


		return result;
	}

	//protected STITree<Object> decompose(STITree<Object>[] U1, Set<String>[] U1_leaf_names, STITree<Object> u2, Tree backbone_tree, Hashtable<STITree<Object>,STITree<Object>> V) {
	public STITree<Object> decompose(LinkedList<STITree<Object>> U1, LinkedList<Set<String>> U1_leaf_names, STITree<Object> u2, Tree backbone_tree, Hashtable<STITree<Object>, STITree<Object>> V) {
		// is there a tree in U1 that contains u2?
		Set<String> u2_leaf_names = getLeafNames(u2);
		for(int i = 0; i < U1.size(); i++) {
			Set<String> s = new HashSet<String>(u2_leaf_names);
			s.removeAll(U1_leaf_names.get(i));

			// if(L(u_2) \subseteq L(u_1) ...
			if(s.isEmpty()) {
				// U' = U' \cup {u_2}
				// V.put(u2, U1.get(i));

				// B(u_2) = T
				_backbone_map.put(u2, backbone_tree);

				/* It can happen that two or more u2 with L(u2) \subset L(u1), but they are not 'independent'.
				 * Clade 'independence' here means that they share some edges in L(u1), so that when there are
				 * SPRs to move those clades, those shared edges make the solution incorrect. The quick fix
				 * is to break u1 into smaller clades such that no two u2 share any edge in one u1. The other
				 * way is to break u2 into smaller piecies so that they do not share edges with each other.
				 *
				 * Questions: Does this way complete?
				 */

				STITree<Object> old_u1 = U1.get(i);
				SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(old_u1);
				Set<STINode<Object>> u2_leaves = getLeaves(old_u1, u2_leaf_names);
				TNode u2lca = lcaFinder.getLCA(u2_leaves);

				if (u2lca == old_u1.getRoot()) {
					STITree<Object> new_u11 = new STITree<Object>(old_u1.getRoot(), true);
					new_u11.constrainByLeaves(u2_leaf_names);
					STITree<Object>[] newTrees = createNonMASTSubtrees(old_u1, new_u11);

					// Remove old u1.
					U1.remove(i);
					U1_leaf_names.remove(i);

					// Add new u1.
					U1.add(new_u11);
					U1_leaf_names.add(u2_leaf_names);

					// Add new sub trees decomposed from u1.
					for (int j = 0; j < newTrees.length; j++) {
						U1.add(newTrees[j]);
						U1_leaf_names.add(getLeafNames(newTrees[j]));
					}

					V.put(u2, new_u11);
				}
				else {
					STITree<Object> new_u11 = new STITree<Object>(u2lca);
					U1.add(new_u11);
					U1_leaf_names.add(getLeafNames(new_u11));

					U1.remove(i);	// Remove old u1.
					U1_leaf_names.remove(i);
					TMutableNode parent = (TMutableNode) u2lca.getParent();
					parent.removeChild((TMutableNode) u2lca, false);

					if (parent.getChildCount() == 1) {	// Binary nodes.
						TMutableNode grandParent = parent.getParent();
						if (grandParent != null) {
							grandParent.removeChild(parent, true);
						}
						else {
							old_u1 = new STITree<Object>(parent.getChildren().iterator().next());
						}
					}
					U1.add(old_u1);	// Old u1 that has removed the subtree containing
					U1_leaf_names.add(getLeafNames(old_u1));

					V.put(u2, new_u11);
				}

//				if (u2lca != old_u1.getRoot()) {	// u2 is a subtree in u1. Remove it from u1.
//					STITree<Object> new_u11 = new STITree<Object>(u2lca);
//					U1.add(new_u11);
//					U1_leaf_names.add(getLeafNames(new_u11));
//
//					U1.remove(i);	// Remove old u1.
//					U1_leaf_names.remove(i);
//					TMutableNode parent = (TMutableNode) u2lca.getParent();
//					parent.removeChild((TMutableNode) u2lca, false);
//
//					if (parent.getChildCount() == 1) {	// Binary nodes.
//						TMutableNode grandParent = parent.getParent();
//						if (grandParent != null) {
//							grandParent.removeChild(parent, true);
//						}
//						else {
//							old_u1 = new STITree<Object>(parent.getChildren().iterator().next());
//						}
//					}
//					U1.add(old_u1);	// Old u1 that has removed the subtree containing
//					U1_leaf_names.add(getLeafNames(old_u1));
//
//					V.put(u2, new_u11);
//				}
//				else {	// Remove all trees not in u2.
//					STITree<Object> new_u11 = new STITree<Object>(old_u1.getRoot(), true);
//					new_u11.constrainByLeaves(u2_leaf_names);
//					STITree<Object>[] newTrees = createNonMASTSubtrees(old_u1, new_u11);
//
//					// Remove old u1.
//					U1.remove(i);
//					U1_leaf_names.remove(i);
//
//					// Add new u1.
//					U1.add(new_u11);
//					U1_leaf_names.add(u2_leaf_names);
//
//					// Add new sub trees decomposed from u1.
//					for (int j = 0; j < newTrees.length; j++) {
//						U1.add(newTrees[j]);
//						U1_leaf_names.add(getLeafNames(newTrees[j]));
//					}
//
//					V.put(u2, new_u11);
//				}

				return u2;
			}
		}

		// otherwise...
		STITree<Object> u1 = null;
		Set<String> u1_leaf_names = null;
		for(int i = 0; i < U1.size(); i++) {
			Set<String> s = getLeafNames(U1.get(i));

			if(shareRoot(u2,s)) {
				u1 = U1.get(i);
				u1_leaf_names = U1_leaf_names.get(i);
				break;
			}
		}

		// if the u2 shares the same root as u2|L(u1) for some u1 \in U1
		if(u1 != null) {
			// t = u_2|L(u_1)
			STITree<Object> t = new STITree<Object>(u2.getRoot(), true);
			t.constrainByLeaves(u1_leaf_names);

			// U' = U' \cup {t}
			V.put(t, u1);

			// B(t) = T
			_backbone_map.put(t, backbone_tree);

			// Foreach x \in X ...
			STITree<Object>[] X = (STITree<Object>[]) createNonMASTSubtrees(u2, t);
			for (int i = 0; i < X.length; i++) {
				decompose(U1, U1_leaf_names, X[i], t, V);
			}

			return t;
		}
		else {
			STINode u2_root = u2.getRoot();
			Iterator<STINode<Object>> i = u2_root.getChildren().iterator();

			// x = Decompose(U1, T_{c1}, T, U')
			STINode<Object> c1 = i.next();

			// Store c1 that can be replaced by c2...ck.

			//_primary_nodes.add(_org_gene_tree.getNode(c1.getID()));
			//_decomp_nodes.add(_org_gene_tree.getNode(u2_root.getID()));


			STITree<Object> child_tree1 = new STITree<Object>(c1, true);
			STITree<Object> x = decompose(U1, U1_leaf_names, child_tree1, backbone_tree, V);

			// For i = 2 to k ...
			while (i.hasNext()) {
				STINode<Object> child = i.next();
				STITree<Object> child_tree = new STITree<Object>(child, true);

				decompose(U1, U1_leaf_names, child_tree, x, V);
			}

			return x;
		}
	}

	private boolean shareRoot(STITree<Object> t1, Set<String> leaf_names) {
		boolean share_root = false;
		TMutableNode prev_child = null;

		Set<String> t1leaves = getLeafNames(t1);
		if (!t1leaves.containsAll(leaf_names)) {
			return false;
		}

		// determine where each leaf is ... or at least until we know
		// there's one on each side of the root
		Iterator<String> i = leaf_names.iterator();
		while(i.hasNext() && !share_root) {
			String name = (String) i.next();

			TMutableNode child = t1.getNode(name);

			// if the leaf doesn't even exist in t1, don't try it
			if(child == null) {
				continue;
			}

			// find the ancestor of this child that is directly beneath the root
			TMutableNode node = child;
			while(node.getParent() != t1.getRoot()) {
				node = node.getParent();
			}

			if(prev_child == null) {
				prev_child = node;
			} else if(prev_child != node) {
				share_root = true;
			}
		}

		return share_root;
	}

	private void addSingleHGT(STITree<Object> species_tree, Tree gene_tree, STITree<Object> u2, STITree<Object>[] U2, Tree mast_tree) {

		// Q = L(u_2) \cup L(B(u_2))
		Set<String> Q = new HashSet<String>();
		Set<String> u2_leaf_names = getLeafNames(u2);
		Set<String> bb_leaf_names = getLeafNames(_backbone_map.get(u2));
		Q.addAll(u2_leaf_names);
		Q.addAll(bb_leaf_names);

		// T'' = GT|Q
		STITree<Object> tpp = new STITree<Object>(gene_tree.getRoot(), true);
		tpp.constrainByLeaves(Q);

		// p = lca_{T''}(L(u_2))
		SchieberVishkinLCA lca = new SchieberVishkinLCA(tpp);
		TMutableNode p;
		Set<STINode<Object>> u2_leaves_in_tpp = getLeaves(tpp, u2_leaf_names);

		p = (STINode<Object>) lca.getLCA(u2_leaves_in_tpp);

		// tq = lca_{ST}(L(u_2))
		Set<STINode<Object>> u2_leaves_in_st = getLeaves(species_tree, u2_leaf_names);
		SchieberVishkinLCA st_lca = new SchieberVishkinLCA(species_tree);
		STINode<Object> tq = (STINode<Object>) st_lca.getLCA(u2_leaves_in_st);

		// te = inedge_{ST}(tq)
		TNode te = _org_species_tree.getNode(tq.getID());
		TNode gene_te = _org_gene_tree.getNode(u2.getRoot().getID());	// te in gt.

		// if p is a child of r(T'') and |L(B(u_2))| > 1 ...
		if(p.getParent() == tpp.getRoot() && bb_leaf_names.size() > 1) {
			TNode sq = st_lca.getLCA(getLeaves(species_tree, bb_leaf_names));
			sq = _org_species_tree.getNode(sq.getID());

			HgtEvent event = new HgtEvent(sq, te);

			if (_primary_nodes.contains(gene_te)) {
				event.setPrimaryNode(gene_te);
			}
			else if (getDecompNode(gene_te) != null) {
				event.addDecompNode(gene_te);
			}

			_hgt_events.add(event);

		} else {

			// O = union (L(T_p'))
			Set<STINode<Object>> O = new HashSet<STINode<Object>>();
			Iterator i = p.getParent().getChildren().iterator();
			while(i.hasNext()) {
				TMutableNode child = (TMutableNode) i.next();

				if(child != p) {
					O.addAll(getLeaves(species_tree, getLeafNames(new STITree(child, true))));
				}
			}

			// sq = lca_{ST}(O)
			TNode sq = st_lca.getLCA(O);

			// se = inedge_{ST}(sq)
			TNode se = _org_species_tree.getNode(sq.getID());

			HgtEvent event = new HgtEvent(se, te);
			if (_primary_nodes.contains(gene_te)) {
				event.setPrimaryNode(gene_te);
			}
			else if (getDecompNode(gene_te) != null) {
				event.addDecompNode(gene_te);
			}

			_hgt_events.add(event);
		}
	}

	/**
	 * Retrieves the set of nodes from tree <code>tree</code> with names contained in set
	 * <code>names</code>.
	 */
	protected Set<STINode<Object>> getLeaves(STITree<Object> tree, Set<String> names) {

		Set<STINode<Object>> s = new HashSet<STINode<Object>>();

		for(String name : names) {
			s.add(tree.getNode(name));
		}

		return s;
	}

	public String printEvents()
	{
		String str = new String();
		boolean hasBad = false;

		int count = 0;
		for (HgtEvent event : _hgt_events) {
			if (!event.isBad()) {
				count++;
				str += (event.getSourceEdge().getName() + " -> " + event.getDestEdge().getName());
				if (event.isViolated()) {
					str += " [time violation?]\n";
				}
				else {
					str += "\n";
				}
			}
			else
				hasBad = true;
		}

		if (hasBad) {
			str += "[refine nodes: ";
			for (HgtEvent event : _hgt_events) {
				if (event.isBad()) {
					str += event.getSourceEdge().getName() + ", ";
				}
			}

			str = str.substring(0, str.length() - 2) + "]\n";
		}
		str = "Number of events: " + count + "\n" + str;
		return str;
	}
}

/**
 * This class is for computing multiple Riata solutions.
 */
class MultipleRiataHgt extends SingleRiataHgt {
	/**
	 * Default constructor: Initializes the class data members.
	 */
	protected MultipleRiataHgt()
	{
		super();
		_hgt_result_tree = null;
		_min_hgt_list = null;
		_min_solution_size = Integer.MAX_VALUE;
	}

	/**
	 * Count the number of HGT scenarios with minimum number of HGT events.
	 *
	 * @return: Number of minimum solutions.
	 */
	public int countMinScenarios()
	{
		if (_min_solution_size == Integer.MAX_VALUE) {
			getMinimumSolutionSize();
		}
		if (_num_min_solution == Integer.MAX_VALUE) {
			getMinHgtScenarios();
		}

		return _min_hgt_list.size();
	}

	/**
	 * This function returns all minimal HGT scenarios.
	 *
	 *  @return A list of minimal HGT scenarios.
	 */
	public LinkedList<HgtScenario> getMinHgtScenarios()
	{
		if (_min_solution_size == Integer.MAX_VALUE) {
			getMinimumSolutionSize();
		}

		if (_min_solution_size == 0) {
			_num_min_solution = 0;
			return new LinkedList<HgtScenario>();	// Return empty list.
		}
		else {
			LinkedList<HgtScenario> scenarioList = getAllHgtScenarios(_hgt_result_tree.getRoot());
			Iterator<HgtScenario> it = scenarioList.iterator();
			while (it.hasNext()) {
				HgtScenario scenario = it.next();

				if (scenario.getEvents().size() - scenario.countBadEvents() != _min_solution_size) {
					// Not the minimal scenario, delete it.
					it.remove();
				}
			}

			_num_min_solution = scenarioList.size();
			return scenarioList;
		}
	}

	/**
	 * This function gets all HGT events from the _hgtResultTree.
	 *
	 * @param node: Indicates the start node in the _hgtResultTree. To get all HGT events
	 * for, call this function with node as the root of _hgtResultTree.
	 *
	 * @return: A linked list of ExHgtScenario.
	 */
	private LinkedList<HgtScenario> getAllHgtScenarios(STINode<ExTNodeData> node)
	{
		// Single-leaf case.
		if (node.isLeaf()) {
			LinkedList<HgtScenario> scenarioList = new LinkedList<HgtScenario>();
			HgtScenario scenario = new HgtScenario();

			if (!node.getData()._mast_node) {	// Event node
				scenario.getEvents().add(node.getData()._event);
			}
			scenarioList.add(scenario);

			return scenarioList;
		}

		// Internal node--there're two cases: MAST node and ordinary node.
		if (node.getData()._mast_node) {	// MAST node.
			// Get scenario lists from this node's children.
			Set<LinkedList<HgtScenario>> s = new HashSet<LinkedList<HgtScenario>>();

			for (STINode<ExTNodeData> child : node.getChildren()) {
				LinkedList<HgtScenario> temp = getAllHgtScenarios(child);

				if (temp != null) {
					s.add(temp);
				}
			}

			// Build the final list of HGT events from the set temp.
			return getAllHgtScenariosHelper(s);
		}
		else {	// Event node
			// Get events from the node's children.
			LinkedList<HgtScenario> scenarioList = new LinkedList<HgtScenario>();

			for (STINode<ExTNodeData> child : node.getChildren()) {
				LinkedList<HgtScenario> temp = getAllHgtScenarios(child);
				if (temp != null) {
					scenarioList.addAll(temp);
				}
			}

			// Add the event at this node to all scenario's in the list.
			HgtEvent event = node.getData()._event;
			if (event != null) {
				for (HgtScenario scenario : scenarioList) {
					scenario.getEvents().add(event);
				}
			}

			return scenarioList;
		}
	}

	/**
	 * This function is used by function <code>getAllHgtScenarios</code> to find all HGT scenarios.
	 *
	 * @return: A list of HGT scenarios.
	 */
	private LinkedList<HgtScenario> getAllHgtScenariosHelper(Set<LinkedList<HgtScenario>> s)
	{
		LinkedList<HgtScenario> result = new LinkedList<HgtScenario>();
		LinkedList<HgtScenario> temp = new LinkedList<HgtScenario>();

		for (LinkedList<HgtScenario> list : s) {
			temp.clear();
			temp.addAll(result);	// Store the result for augmenting.
			result.clear();

			if (temp.isEmpty()) {
				result.addAll(list);
			}
			else {
				Iterator tempIt = temp.iterator();
				while (tempIt.hasNext()) {
					HgtScenario tempScenario = (HgtScenario) tempIt.next();
					Iterator listIt = list.iterator();
					while (listIt.hasNext()) {
						HgtScenario scenario = new HgtScenario();
						scenario.getEvents().addAll(((HgtScenario) listIt.next()).getEvents());
						scenario.getEvents().addAll(tempScenario.getEvents());

						result.add(scenario);
					}
				}
			}
		}

		return result;
	}

	/**
	 * This function computes ALL possible HGT scenarios for the two species and gene trees.
	 *
	 * @param speciesTree: The species tree.
	 * @param geneTree: The gene tree.
	 */
	public void computeMultipleHgt(Tree speciesTree, Tree geneTree)
	{
		if(!speciesTree.isRooted() || !geneTree.isRooted()) {
			throw new RuntimeException("Trees must be rooted");
		}

		_hgt_result_tree = new STITree<ExTNodeData>();
		_hgt_result_tree.getRoot().setData(new ExTNodeData(false, 0, null, null));

		_org_gene_tree = geneTree;
		_org_species_tree = speciesTree;

		computeMultipleHgtHelper(speciesTree, geneTree, _hgt_result_tree.getRoot());
	}

	/**
	 * This function computes the size of a minimal solutions, i.e. the number
	 * of HGT events in a minimal solutions. Note that the size of a minimal solutions
	 * does not count toward events like x -> y where x is a parent of y. These events
	 * are called bad events and can be substituted by refinement on x.
	 *
	 * @return The size of a minimal solutions.
	 */
	public int getMinimumSolutionSize()
	{
		_min_solution_size = getMinimumSolutionSize(_hgt_result_tree.getRoot());
		return _min_solution_size;
	}

	/**
	 * This recursive function is called by getMinimumSolutionSize.
	 */
	private int getMinimumSolutionSize(STINode<ExTNodeData> node)
	{
		if (node.isLeaf()) {
			if (node.getData()._mast_node) {	// MAST node--no event here.
				return 0;
			}
			else {	// Ordinary node. There's exactly one event.
				HgtEvent event = node.getData()._event;

				if (event == null || event.isBad()) {
					return 0;	// Do not count bad events.
				}
				else {
					return 1;
				}
			}
		}

		int count;
		if (node.getData()._mast_node) {	// A MAST node.
			// The min # of events is equal to the sum of of events of its children.
			count = 0;

			for (STINode<ExTNodeData> child : node.getChildren()) {
				if (child.getData()._mast_node) {
					System.out.println("A MAST node will never have a child being also a MAST.");
					assert(false);
				}
				count += getMinimumSolutionSize(child);
			}

			return count;
		}
		else {	// An event node.
			// Find the minium solution size among its children.
			count = Integer.MAX_VALUE;

			for (STINode<ExTNodeData> child : node.getChildren()) {
				if (!child.getData()._mast_node) {
					System.out.println("An ordinary node cannot have children as ordinary nodes.");
					assert(false);
				}

				int temp = getMinimumSolutionSize(child);
				if (count > temp) {
					count = temp;
				}
			}

			// Plus one event stored within this event node.
			HgtEvent event = node.getData()._event;

			if (event == null || event.isBad()) {
				return count;	// Do not count bad events.
			}
			else {
				return count + 1;
			}
		}
	}

	private void computeMultipleHgtHelper(Tree speciesTree, Tree geneTree, STINode<ExTNodeData> parent)
	{
		// compute the rooted mast for trees with an outgroup. (T')
		ExMultipleMasts emm = new ExMultipleMasts();
		Tree xst = SingleRiataHgt.addOutgroup(speciesTree);
		Tree xgt = SingleRiataHgt.addOutgroup(geneTree);

		Set<Tree> xmts = emm.computeMultipleMasts(xst, xgt);	//Set<Tree> mastTreeSet = emm.exComputeMultipleMasts(speciesTree, geneTree);
		if (xmts.isEmpty()) {
			return;	// The two trees are identical.
		}

		// Remove the outgroup from MASTs
		Set<Tree> mastTreeSet = new HashSet<Tree>();

		for (Tree xmast : xmts) {
			Tree mast = SingleRiataHgt.removeOutgroup(xmast);
			Set<String> mastleafset = getLeafNames(mast);
			mastTreeSet.add(mast);
		}

		// Compute HGT events.
		STITree<Object> sti_st = new STITree<Object>(speciesTree.getRoot(), true);
		for (Tree mastTree : mastTreeSet) {
			STINode<ExTNodeData> mastNode = parent.createChild();
			mastNode.setData(new ExTNodeData(true, parent.getData()._depth + 1, null, mastTree));


			// If the trees don't overlap at all, there's no HGT information encoded here.
			if(mastTree == null || mastTree.getLeafCount() == 0) {
				return;
			}

			// If the mast is the species tree, there's no difference between the trees
			if(mastTree.getLeafCount() == speciesTree.getLeafCount()) {
				return;
			}

			// U1 = ST - T'
			STITree<Object>[] U1_array = createNonMASTSubtrees(speciesTree, mastTree);
			LinkedList<STITree<Object>> U1 = new LinkedList<STITree<Object>>();
			for (int i = 0; i < U1_array.length; i++) {
				U1.add(U1_array[i]);
			}

			LinkedList<Set<String>> U1_leaf_names = new LinkedList<Set<String>>();
			for(int i = 0; i < U1.size(); i++) {
				U1_leaf_names.add(getLeafNames(U1.get(i)));
			}

			// U2 = GT - T'
			STITree<Object>[] U2 = createNonMASTSubtrees(geneTree, mastTree);

			// Decompose the clades in u2
			// Store decomposition nodes for recovery later.
			List<TNode> decomp_copies = new LinkedList<TNode>();
			List<TNode> primary_copies = new LinkedList<TNode>();
			decomp_copies.addAll(_decomp_nodes);
			primary_copies.addAll(_primary_nodes);

			Hashtable<STITree<Object>, STITree<Object>> V = new Hashtable<STITree<Object>, STITree<Object>>();

			for(int i = 0; i < U2.length; i++) {
				decompose(U1, U1_leaf_names, U2[i], mastTree, V);
			}

			U2 = (STITree<Object>[]) new STITree[V.size()];
			Iterator iter = V.keySet().iterator();
			for(int i = 0; iter.hasNext(); i++) {
				U2[i] = (STITree<Object>) iter.next();
			}

			// while ...
			for (int i = 0; i < U2.length; i++) {
				STITree<Object> u2 = U2[i];
				STITree<Object> u1 = V.get(u2);	// Corresponding u1 s.t. L(u2) is contained in L(u1).
				Set<String> u2_leafset = getLeafNames(u2);
				Set<String> u1_leafset = getLeafNames(u1);

				// DEBUG
				Set<String> temp = new HashSet<String>();
				temp.addAll(u2_leafset);
				temp.removeAll(u1_leafset);
				if (!temp.isEmpty()) {
					System.err.println("u2: " + u2);
					System.err.println("u1: " + u1);
					System.err.println("Leaves in u2 are not contained in u1");
					return;
				}
				// DEBUG

				STITree<Object> u1_copy = new STITree<Object>(u1.getRoot(), true);
				u1_copy.constrainByLeaves(u2_leafset);

				// Store the backbone map so that we can continue to find new events.
				Hashtable oldBackboneMap = new Hashtable(_backbone_map);

				// Compute HGT events recursively.
				STINode<ExTNodeData> normalNode = mastNode.createChild();
				normalNode.setData(new ExTNodeData(false, mastNode.getData()._depth + 1, null, null));

				computeMultipleHgtHelper(u1_copy, u2, normalNode);
				addSingleHgt(sti_st, (STITree<Object>) geneTree, u2, U2, mastTree, normalNode);

				// Restore information here.
				_backbone_map.clear();
				_backbone_map.putAll(oldBackboneMap);
			}

			// Restore information about decomposition nodes.
			_decomp_nodes.clear();
			_decomp_nodes.addAll(decomp_copies);
			_primary_nodes.clear();
			_primary_nodes.addAll(primary_copies);
		}
	}

	private void addSingleHgt(STITree<Object> species_tree, STITree<Object> gene_tree, STITree<Object> u2, STITree<Object>[] U2, Tree mast_tree, STINode<ExTNodeData> eventNode)
	{
		// Q = L(u_2) \cup L(B(u_2))
		Set<String> Q = new HashSet<String>();
		Set<String> u2_leaf_names = getLeafNames(u2);
		Set<String> bb_leaf_names = getLeafNames((STITree) _backbone_map.get(u2));
		Q.addAll(u2_leaf_names);
		Q.addAll(bb_leaf_names);

		// T'' = GT|Q
		STITree tpp = new STITree(gene_tree.getRoot(), true);
		tpp.constrainByLeaves(Q);

		// p = lca_{T''}(L(u_2))
		SchieberVishkinLCA lca = new SchieberVishkinLCA(tpp);
		TNode p;
		Set u2_leaves_in_tpp = getLeaves(tpp, u2_leaf_names);

		p = lca.getLCA(u2_leaves_in_tpp);

		// tq = lca_{ST}(L(u_2))
		Set u2_leaves_in_st = getLeaves(species_tree, u2_leaf_names);
		SchieberVishkinLCA st_lca = new SchieberVishkinLCA(species_tree);
		TNode tq = st_lca.getLCA(u2_leaves_in_st);

		// te = inedge_{ST}(tq)
		TNode te = _org_species_tree.getNode(tq.getID());
		TNode gene_te = _org_gene_tree.getNode(u2.getRoot().getID());	// te in gt.


		// if p is a child of r(T'') and |L(B(u_2))| > 1 ...
		if(p.getParent() == tpp.getRoot() && bb_leaf_names.size() > 1) {
			TNode sq = st_lca.getLCA(getLeaves(species_tree, bb_leaf_names));
			HgtEvent event = new HgtEvent(sq, te);

			if (_primary_nodes.contains(gene_te)) {
				event.setPrimaryNode(gene_te);
			}
			else if (getDecompNode(gene_te) != null) {
				event.addDecompNode(gene_te);
			}

			eventNode.getData().setHgtEvent(event);
		}
		else {
			// O = union (L(T_p'))
			Set O = new HashSet();
			Iterator i = p.getParent().getChildren().iterator();
			while(i.hasNext()) {
				TNode child = (TNode) i.next();

				if(child != p) {
					O.addAll(getLeaves(species_tree, getLeafNames(new STITree(child, true))));
				}
			}

			TNode sq = st_lca.getLCA(O);	// sq = lca_{ST}(O)
			TNode se = _org_species_tree.getNode(sq.getID());

			HgtEvent event = new HgtEvent(se, te);

			if (_primary_nodes.contains(gene_te)) {
				event.setPrimaryNode(gene_te);
			}
			else {
				TNode decomp_node = getDecompNode(gene_te);
				if (decomp_node != null) {
					event.addDecompNode(decomp_node);
				}
			}

			eventNode.getData().setHgtEvent(event);
		}
	}

	/**
	 * This function clear the name of internal nodes. This function is used
	 * to clear names in the gene tree to avoid name conflict.
	 */
	private static void clearNodeNames(MutableTree tree) {
		for(TMutableNode n : tree.getNodes()) {
			if (!n.isLeaf())
				n.setName(TNode.NO_NAME);
		}
	}


	/* Data members. */
	public STITree<ExTNodeData> _hgt_result_tree = null;	// The tree that stores all HGT events.
	public LinkedList<HgtScenario> _min_hgt_list;			// The list of HGT scenarios with minimum number of events. Each element is an ExHgtScenario.
	public int _min_solution_size = Integer.MAX_VALUE;	// The minimum number of events.
	public int _num_min_solution = Integer.MAX_VALUE;	// The number of minimal solutions.

	/* Data structure to store information inside a node in the result tree. */
	private class ExTNodeData {
		/**
		 * This constructor build an instance of ExTNodeData.
		 *
		 * @param is_mast = <code>true</code> if we want to create a node to store
		 * a MAST tree; it is <code>false</code> if this instance is for storing an event.
		 * @param depth: is the depth of this node in the result tree. This value is for
		 * debugging only. It can be safely set to zero.
		 * @param event: If this node is to store an event, then <code>event</code> should not
		 * be null.
		 * @param mast: If this node is to store a MAST, then <code>mast</code> should not be null.
		 */
		ExTNodeData(boolean is_mast, int depth, HgtEvent event, Tree mast)
		{
			_mast_node = is_mast;
			_depth = depth;

			if (_mast_node) {
				_event = null;
				_mast_tree = mast;
			}
			else {
				_event = event;
				_mast_tree = null;
			}
		}

		/**
		 * This constructor build an instance of ExTNodeData.
		 *
		 * @param is_mast = <code>true</code> if we want to create a node to store
		 * a MAST tree; it is <code>false</code> if this instance is for storing an event.
		 * @param depth: is the depth of this node in the result tree. This value is for
		 * debugging only. It can be safely set to zero.
		 * @param src and dest: specify the donor and receipent of an event. These values should
		 * not be null if this node is to store an event.
		 * @param mast: If this node is to store a MAST, then <code>mast</code> should not be null.
		 */
		ExTNodeData(boolean is_mast, int depth, TNode src, TNode dest, STITree mast)
		{
			this(is_mast, depth, new HgtEvent(src, dest), mast);
		}

		/**
		 * This function stores an event in a node of the result tree.
		 *
		 * @param event: The event to be stored in a node.
		 */
		public void setHgtEvent(HgtEvent event)
		{
			_event = event;
			_mast_node = false;
		}

		/**
		 * This function stores an event in a node of the result tree.
		 *
		 * @param src: The donor of an event.
		 * @param dest: The receipent of an event.
		 */
		public void setHgtEvent(TNode src, TNode dest)
		{
			_event = new HgtEvent(src, dest);
			_mast_node = false;
		}

		/**
		 * This function stores a reference to a MAST tree to a node of the
		 * result tree.
		 *
		 * @param mastTree: The reference to the MAST tree to be stored.
		 */
		public void setMastTree(STITree mastTree)
		{
			_mast_tree = mastTree;
			_mast_node = true;
		}

		/**
		 * This method gets a string representing the MAST tree of this MAST node.
		 *
		 * @return: A Newick string representing the MAST stored in this node.
		 */
		public String toMastTree()
		{
			String result = new String();

			if (_mast_node) {
				for (int i = 0; i < _depth; i++) {
					result += "\t";
				}

				if (_mast_tree != null) {
					result += "MAST node: " + _mast_tree.toNewick();
				}
			}

			return result;
		}

		/**
		 * This method gets a string representing the event of this event node.
		 */
		public String toString()
		{
			String result = new String();

			if (!_mast_node) {
				for (int i = 0; i < _depth; i++) {
					result += "\t";
				}

				if (_event != null) {
					result += _event.getSourceEdge().getName() + " -> " + _event.getDestEdge().getName();
				}
			}

			return result;
		}

		/* Data members of ExTNodeData */
		public boolean _mast_node;	// true if this node is to store a MAST.
		public int _depth;			// Stores the depth of this node. It's for debugging purposes.
		public HgtEvent _event;		// The event at an ordinary node (_mast_node = false).
		public Tree _mast_tree;		// The corresponding MAST (not null only if _mast_node = true).
	}
}

