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

package edu.rice.cs.bioinfo.programs.phylonet.algos.consensus;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/30/11
 * Time: 4:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class TreeConsensusCalculator
{
    // constants
	private static String SEP_STR = "###";

	// inner classes
	private class NamesObject {
		public String[] names;
		public double support;

		public NamesObject(String[] names, double support) {
			this.names = names;
			this.support = support;
		}
	}

	// fields
	private Hashtable<String,Double> _tally = new Hashtable<String,Double>();

	private boolean _input_supports = false;
	private boolean _output_supports = false;

	// methods
	/**
	 * If this property is set to <code>true</code>, then edge weights on
	 * input trees are interpretted as existing edge supports.  These are
	 * included in the computation of the overall weights.  If an edge
	 * does not have any weight, then the support is assumed to be 1.
	 * All supports must be between 0 and 1.
	 */
	public void setInputEdgeWeightsAreSupports(boolean supports) {
		_input_supports = supports;
	}

	public boolean getInputEdgeWeightsAreSupports() { return _input_supports; }

	public void setOutputEdgeWeightsAreSupports(boolean supports) {
		_output_supports = supports;
	}

	public boolean getOutputEdgeWeightsAreSupports() { return _output_supports; }

	public MutableTree computeUnrootedConsensus(Tree[] trees, double threshold) {
		return computeUnrootedConsensus(trees, trees.length, threshold);
	}

	public MutableTree computeUnrootedConsensus(Tree[] trees, int num_trees, double threshold) {

		// select the 'root' leaf
		String root_leaf_name = null;
		for(TNode n : trees[0].getNodes()) {
			if(n.isLeaf()) {
				root_leaf_name = n.getName();
				break;
			}
		}

		// reroot the trees at the same leaf
		HashSet<STITree<Object>> ntrees = new HashSet<STITree<Object>>();
		int idx = 0;
		for(Tree t : trees) {

			if(idx == num_trees) {
				break;
			}

			STITree<Object> nt = new STITree<Object>(t.getRoot(), true);
			nt.getNode(root_leaf_name).makeRoot();
			ntrees.add(nt);

			idx++;
		}

		return inner_computeUnrootedConsensus(ntrees, threshold, root_leaf_name);
	}

	public MutableTree computeUnrootedConsensus(Set<? extends Tree> trees, double threshold) {

		// select the 'root' leaf
		String root_leaf_name = null;
		for(TNode n : trees.iterator().next().getNodes()) {
			if(n.isLeaf()) {
				root_leaf_name = n.getName();
				break;
			}
		}

		// reroot the trees at the same leaf
		HashSet<STITree<Object>> ntrees = new HashSet<STITree<Object>>();
		for(Tree t : trees) {
			STITree<Object> nt = new STITree<Object>(t.getRoot(), true);
			nt.getNode(root_leaf_name).makeRoot();
			ntrees.add(nt);
		}

		return inner_computeUnrootedConsensus(ntrees, threshold, root_leaf_name);
	}

	protected MutableTree inner_computeUnrootedConsensus(Set<STITree<Object>> ntrees, double threshold, String root_leaf_name) {

		// perform the rooted consensus
		MutableTree result = computeRootedConsensus(ntrees, threshold);

		// correct for the loss of a leaf by adding the child back in
		// this happens because the consensus algorithm automatically removes
		// extraneous inner nodes that have only one child.
		TMutableNode n = result.getRoot().createChild(root_leaf_name);

		if(_output_supports) {
			n.setParentDistance(1.0);
		}

		// done!
		return result;
	}

	/**
	 * This method computes the consensus of a set of trees.  The edges
	 * retained are those that occur in at least <code>threshold</code> of
	 * the trees.
	 *
	 * @param tree_set is a set of trees
	 * @param threshold is a number between 0 and 1.
	 */
	public MutableTree computeRootedConsensus(Set<? extends Tree> tree_set, double threshold) {

		// make sure that all binary nodes are removed (to avoid double counting edges)
		Set<STITree<Object>> trees = new HashSet<STITree<Object>>();
		for(Tree t : tree_set) {
			STITree<Object> nt = new STITree<Object>(t);
			Trees.removeBinaryNodes(nt);
			trees.add(nt);
		}

		return inner_computeRootedConsensus(trees, threshold);
	}

	/**
	 * Compute the consensus for all the trees in the array.
	 */
	public MutableTree computeRootedConsensus(Tree[] tree_array, double threshold) {
		return computeRootedConsensus(tree_array, tree_array.length, threshold);
	}

	/**
	 * Compute the consensus for only a subset of the trees in the array.
	 */
	public MutableTree computeRootedConsensus(Tree[] tree_array, int num_trees, double threshold) {

		// make sure that all binary nodes are removed (to avoid double counting edges)
		Set<STITree<Object>> trees = new HashSet<STITree<Object>>();
		int idx = 0;
		for(Tree t : tree_array) {

			if(idx == num_trees) {
				break;
			}

			STITree<Object> nt = new STITree<Object>(t);
			Trees.removeBinaryNodes(nt);
			trees.add(nt);

			idx++;
		}

		return inner_computeRootedConsensus(trees, threshold);
	}

	/**
	 * This method does the actual work computing the consensus
	 *
	 * @param trees
	 * @param threshold
	 * @return
	 */
	protected MutableTree inner_computeRootedConsensus(Set<STITree<Object>> trees, double threshold) {

		int idx;

		// Make sure that leaf sets are identical
		Iterator<STITree<Object>> i = trees.iterator();
		MutableTree t1 = i.next();

		while(i.hasNext()) {
			MutableTree t2 = i.next();

			if(!Trees.leafSetsAgree(t1, t2)) {
				throw new RuntimeException("Leafsets do not agree in trees");
			}

			if(!t1.isRooted() || !t2.isRooted()) {
				throw new RuntimeException("Trees must be rooted");
			}
		}

		// count the number of times each edge appears
		_tally.clear();

		for(Tree t : trees) {
			tallyEdges(t);
		}

		// select the edges that pass the threshold
		Iterator<Map.Entry<String,Double>> it = _tally.entrySet().iterator();
		while(it.hasNext()) {
			Map.Entry<String,Double> e = it.next();

			double support = ((Double) e.getValue()).doubleValue();
			double percent = (double) support / (double) trees.size();

			if(percent < threshold) {
				it.remove();
			} else {
				e.setValue(new Double(percent));
			}
		}

		// make a list of only those edges that pass the threshold
		NamesObject[] names_array = new NamesObject[_tally.size()];
		idx = 0;
		for(Map.Entry<String,Double> e : _tally.entrySet()) {
			String name = (String) e.getKey();
			double support = ((Double) e.getValue()).doubleValue();

			names_array[idx++] = new NamesObject(edge2LeafNames(name), support);
		}

		// order the edges by number of leaves each has under it
		Arrays.sort(names_array, new Comparator<Object>() {

			public int compare(Object o1, Object o2) {
				NamesObject nlist1 = (NamesObject) o1;
				NamesObject nlist2 = (NamesObject) o2;

				if(nlist1.names.length < nlist2.names.length) {
					return -1;
				} else if(nlist1.names.length == nlist2.names.length) {
					return 0;
				} else {
					return 1;
				}
			}
		});

		// build the new tree from these edges starting with the smallest ones and going up
		STITree<Object> result_tree = new STITree<Object>();
		Hashtable<String,STINode<Object>> leaf2node = new Hashtable<String,STINode<Object>>();

		for(idx = 0; idx < (names_array.length - 1); idx++) {
			// create a new node for this edge
			STINode<Object> node = result_tree.getRoot().createChild();

			if(_output_supports) {
				node.setParentDistance(names_array[idx].support);
			}

			String[] names = names_array[idx].names;

			// the leaf case
			if(names.length == 1) {
				node.setName(names[0]);
				leaf2node.put(names[0], node);
			} else {
				// add all the nodes beneath it
				for(String name : names) {
					TMutableNode child = leaf2node.get(name);
					node.adoptChild(child);
					leaf2node.put(name, node);
				}
			}
		}

		// remove binary nodes
		Trees.removeBinaryNodes(result_tree);

		// done!
		return result_tree;
	}

	private void tallyEdges(Tree t) {
		tallyEdges(t.getRoot());
	}

	private Set<TNode> tallyEdges(TNode node) {

		Set<TNode> s = null;

		if(node.isLeaf()) {
			s = new HashSet<TNode>();
			s.add(node);
		} else {
			Iterator i = node.getChildren().iterator();
			s = tallyEdges((TMutableNode) i.next());

			while(i.hasNext()) {
				s.addAll(tallyEdges((TMutableNode) i.next()));
			}
		}

		// update the count of the edge
		String name = leafSet2EdgeName(s);
		if(!_tally.containsKey(name)) {
			if(!_input_supports || node.getParentDistance() == TNode.NO_DISTANCE) {
				_tally.put(name, 1.0);
			} else {
				_tally.put(name, node.getParentDistance());
			}
		} else {
			if(!_input_supports || node.getParentDistance() == TNode.NO_DISTANCE) {
				_tally.put(name, new Double(((Double) _tally.get(name)).doubleValue() + 1.0));
			} else {
				_tally.put(name, new Double(((Double) _tally.get(name)).doubleValue() + node.getParentDistance()));
			}
		}

		return s;
	}

	private String leafSet2EdgeName(Set leaf_set) {
		String name = "";

		// sort the leaves by their name
		Object[] array = leaf_set.toArray();

		Arrays.sort(array, new Comparator<Object>() {

			public int compare(Object o1, Object o2) {

				TNode n1 = (TNode) o1;
				TNode n2 = (TNode) o2;

				return n1.getName().compareTo(n2.getName());
			}

		});

		// make one string
		name = ((TNode) array[0]).getName();
		for(int i = 1; i < array.length; i++) {
			name = name + SEP_STR + ((TNode) array[i]).getName();
		}

		return name;
	}

	private String[] edge2LeafNames(String str) {
		return str.split(SEP_STR);
	}
}
