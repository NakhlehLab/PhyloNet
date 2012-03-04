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

package edu.rice.cs.bioinfo.programs.phylonet.algos.mast;

import edu.rice.cs.bioinfo.programs.phylonet.algos.bipartitematching.HungarianBipartiteMatcher;
import edu.rice.cs.bioinfo.programs.phylonet.structs.BitVector;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/* This class implements the algorithm for solving the Maximum Agreement Subtree Problem (MAST) described
* in the paper by Mike Steel and Tandy Warnow entitled <B>Kaikoura tree theorems: Computing the maximum
* agreement subtree</B>.
*
* @author Derek Ruths
*/
public class SteelWarnowMAST {

	// inner classes
	/**
	 * This class represents a pairing to two e-subtrees: one from one tree, one from another tree.
	 * This class exists primarily for making it easy to store the results of subtree MAST calculations.
	 * It can be used as a key to find the MAST of these two e-subtrees.  For this to work, the
	 * convention of <tree1,tree2> must be maintained - that is, the order in which the trees are
	 * specified can't change.
	 *
	 * @author Derek Ruths
	 *
	 */
	protected static class ESubtreePair {

		BitVector leaf_vector1;
		BitVector leaf_vector2;

		public ESubtreePair() {}

		public ESubtreePair(BitVector lv1, BitVector lv2) {
			set(lv1,lv2);
		}

		/**
		 * Constructor to build an esubtree pair from two esubtrees e1 and e2.
		 *
		 * @param e1
		 * @param e2
		 *
		 * NOTE: Added by Cuong Than on Mar. 18, 06.
		 */
		public ESubtreePair(ESubtree e1, ESubtree e2)
		{
			this(e1.leaf_vector, e2.leaf_vector);
		}

		public String toString() {
			return "(" + leaf_vector1.toString() + "," + leaf_vector2.toString() + ")";
		}

		public void set(ESubtree s1, ESubtree s2) {
			this.leaf_vector1 = s1.leaf_vector;
			this.leaf_vector2 = s2.leaf_vector;
		}

		public void set(BitVector lv1, BitVector lv2) {
			this.leaf_vector1 = lv1;
			this.leaf_vector2 = lv2;
		}

		/**
		 * This is a simple hash code for this object.  It might be worth looking
		 * at this to see if there is a better hashCode that would decrease
		 * hash collisions.
		 */
		public int hashCode() {
			return leaf_vector1.hashCode() + leaf_vector2.hashCode();
		}

		public boolean equals(Object obj) {
			ESubtreePair p = (ESubtreePair) obj;

			return (p.leaf_vector1.equals(leaf_vector1) && p.leaf_vector2.equals(leaf_vector2));
		}
	}

	/**
	 * This class represents a subtree of tree formed by cutting one of its edges and rooting it at a
	 * node on either side of the cut.  Because e-subtrees will only be comapred with other e-subtrees
	 * created from the same tree, only the leaf-set of the e-subtree is required to distinguish it
	 * from others.  As a result, the only data stored inside the e-subtree for identification purposes
	 * is the leaf-vector.
	 *
	 * @author Derek Ruths
	 */
	protected static class ESubtree {
		/**
		 * This is all that is required to distinguish this e-subtree from other e-subtrees formed
		 * from the same tree.
		 */
		public BitVector leaf_vector;

		/**
		 * These are the (unique) immediate subtrees of the root of
		 * this e-subtree
		 */
		public ESubtree[] subtrees;

		/**
		 * This flag indicates whether this e-subtree is involved in the final computation of a rooted mast of
		 * the tree from which it was derived.
		 */
		public boolean is_root;

		public int num_refs = 0;

		/**
		 * When the e-subtree is formed by removing an edge, another e-subtree is necessarily formed (the tree
		 * has been divided in two parts).  The complement is the tree on the other side of the cut.
		 */
		public ESubtree complement;

		public ESubtree(SteelWarnowMAST swmast, BitVector leaves, ESubtree[] subtrees) {
			this.subtrees = subtrees;
			leaf_vector = leaves;
		}

		public int getRefCount() { return num_refs; }

		public void incrementRefs() {
			num_refs++;
		}

		public void decrementRefs() {

			if(num_refs == 0) {
				throw new RuntimeException("Cannot decrement zero refs");
			}

			num_refs--;
		}

		public boolean canDeleteRMASTInfo() { return (num_refs == 0); }
	}

	// fields
	protected BitVector EMPTY_SET = null;

	protected Hashtable<String,Integer> _lname2num;
	protected String[] _lnum2name;

	protected int _num_leaves;

	protected ESubtree[] _esubtrees1;
	protected ESubtree[] _esubtrees2;

	protected Hashtable<BitVector,ESubtree> _leafset2subtree;

	protected Hashtable<ESubtreePair,BitVector> _mast_results;

	// constructors

	// methods
	/**
	 * This method computes the rooted mast for a set of trees.
	 */
	public Tree computeRMAST(Iterable<? extends Tree> trees) {

		Tree mast = null;

		for(Tree t : trees) {
			if(mast == null) {
				mast = t;
			} else {
				mast = computeRMAST(mast, t);
			}
		}

		return mast;
	}

	/**
	 * This method computes the unrooted mast for a set of trees.
	 */
	public Tree computeUMAST(Iterable<? extends Tree> trees) {

		Tree mast = null;

		for(Tree t : trees) {
			if(mast == null) {
				mast = t;
			} else {
				mast = computeUMAST(mast, t);
			}
		}

		return mast;
	}

	/**
	 * This method computes the rooted maximum agreement subtree for two trees.
	 *
	 * @param t1 is the first tree
	 * @param t2 is the second tree
	 *
	 * @return the maximum agreement subtree.  The tree is not guaranteed to be of the
	 * same type as the input trees.
	 *
	 * @throws RuntimeException if the trees are not both rooted.
	 */
	public Tree computeRMAST(Tree tr1, Tree tr2) {
		// Add an outgroup to tr1 and tr2.
		Tree mast = inner_computeRMAST(tr1, tr2, true);

		return mast;
	}

	/**
	 * Computes the max. agreement subtree for two rooted trees.  This method is
	 * intended for use by implementations extending this base implementation
	 * of the algorithm.  In particular, this function exposes the ability to
	 * keep the method from cleaning up the <code>_mast_result</code> data structure
	 * which contains all submasts computed during the computation.
	 *
	 * @param tr1
	 * @param tr2
	 * @param cleanup if true, intermediate results in <code>_mast_result</code> will
	 * be discarded when they are no longer needed.
	 */
	protected Tree inner_computeRMAST(Tree tr1, Tree tr2, boolean cleanup) {

		/*
		if(!tr1.isRooted() || !tr2.isRooted()) {
			throw new RuntimeException("computeRMAST: Unrooted trees specified");
		}
		*/

		STITree<Object> t1 = new STITree<Object>(tr1);
		STITree<Object> t2 = new STITree<Object>(tr2);

		makeLeafSetsEqual(t1,t2);

		// enforce size
		if(t1.getNodeCount() == 0 || t2.getNodeCount() == 0) {
			// return an empty tree
			return t1;
		}

		/*
		if(!leafSetsEquals(t1,t2)) {
			throw new RuntimeException("computeRMAST: Trees must have identical leafsets");
		}
		*/

		// handle single node tree cases
		if(t1.getNodeCount() == 1) {
			STITree<Object> t = new STITree<Object>(t2.getRoot(), true);

			if(t.getNode(t1.getRoot().getName()) != null) {
				return new STITree<Object>(t1.getRoot(), true);
			} else {
				// create an empty tree
				STITree<Object> tree = new STITree<Object>(true);
				tree.getRoot().removeNode();
				return tree;
			}
		} else if(t2.getNodeCount() == 1) {
			STITree<Object> t = new STITree<Object>(t1.getRoot(), true);

			if(t.getNode(t2.getRoot().getName()) != null) {
				return new STITree<Object>(t2.getRoot(), true);
			} else {
				// create an empty tree
				STITree<Object> tree = new STITree<Object>(true);
				tree.getRoot().removeNode();
				return tree;
			}
		}

		// create unrooted versions of the tree
		STITree<ESubtree> tree1 = new STITree<ESubtree>(t1.getRoot(), false);
		STITree<ESubtree> tree2 = new STITree<ESubtree>(t2.getRoot(), false);

		computeSubMASTs(tree1, tree2, true, cleanup);

		// do the matching between subtrees of the two roots
		ESubtree[] strees1 = new ESubtree[tree1.getRoot().getChildCount()];
		ESubtree[] strees2 = new ESubtree[tree2.getRoot().getChildCount()];

		int i = 0;
		for(STINode<ESubtree> n : tree1.getRoot().getChildren()) {
			strees1[i++] = n.getData();
		}

		i = 0;
		for(STINode<ESubtree> n : tree2.getRoot().getChildren()) {
			strees2[i++] = n.getData();
		}

		BitVector result = calculateMaxWeightedBipartiteMatch(strees1, strees2, _mast_results);

		LinkedList<String> leaf_names = new LinkedList<String>();
		int pos = 0;
		for(boolean bit : result) {
			if(bit) {
				leaf_names.add(_lnum2name[pos]);
			}

			pos++;
		}

		STITree<Object> tresult = new STITree<Object>(tree1.getRoot(), true);
		tresult.constrainByLeaves(leaf_names);

		return tresult;
	}

	/**
	 * This method computes the unrooted maximum agreement subtree for two trees.
	 *
	 * @param t1 is the first tree
	 * @param t2 is the second tree
	 *
	 * @return the maximum agreement subtree.  The tree is not guaranteed to be of the
	 * same type as the input trees.
	 *
	 * @throws RuntimeException if the trees are not both unrooted.
	 */
	public Tree computeUMAST(Tree t1, Tree t2) {

		/*
		if(t1.isRooted() || t2.isRooted()) {
			throw new RuntimeException("computeUMAST: Rooted trees specified");
		}
		*/

		STITree<ESubtree> tree1 = new STITree<ESubtree>(t1.getRoot(), false); //t1.isRooted());
		STITree<ESubtree> tree2 = new STITree<ESubtree>(t2.getRoot(), false); //t2.isRooted());

		// prune the trees to make them equal
		makeLeafSetsEqual(tree1,tree2);

		// check for empty trees
		if(tree1.getNodeCount() == 0 || tree2.getNodeCount() == 0) {
			// return an empty tree
			return tree1;
		}

		/*
		if(!leafSetsEquals(tree1,tree2)) {
			throw new RuntimeException("computeRMAST: Trees must have identical leafsets");
		}
		*/
		// handle single node tree cases
		if(t1.getNodeCount() == 1) {
			STITree<Object> t = new STITree<Object>(t2.getRoot(), true);

			if(t.getNode(t1.getRoot().getName()) != null) {
				return new STITree<Object>(t1.getRoot(), true);
			} else {
				// create an empty tree
				STITree<Object> tree = new STITree<Object>(true);
				tree.getRoot().removeNode();
				return tree;
			}
		} else if(t2.getNodeCount() == 1) {
			STITree<Object> t = new STITree<Object>(t1.getRoot(), true);

			if(t.getNode(t2.getRoot().getName()) != null) {
				return new STITree<Object>(t2.getRoot(), true);
			} else {
				// create an empty tree
				STITree<Object> tree = new STITree<Object>(true);
				tree.getRoot().removeNode();
				return tree;
			}
		}

		// otherwise, do the whole calculation
		computeSubMASTs(tree1, tree2, false, true);

		// compute MAX(p,q)
		ESubtree best_p1 = null;
		ESubtree best_q1 = null;
		int max_size = -1;

		ESubtreePair name = new ESubtreePair();
		ESubtreePair cname = new ESubtreePair();
		BitVector s1;
		BitVector s2;

		for(int i = 1; i < _esubtrees1.length; i++) {

			ESubtree p1 = _esubtrees1[i];

			for(int j = 0; j < _esubtrees2.length; j++) {

				// check p1 with q1
				ESubtree q1 = _esubtrees2[j];

				name.set(p1.leaf_vector,q1.leaf_vector);
				cname.set(p1.complement.leaf_vector,q1.complement.leaf_vector);

				s1 = _mast_results.get(name);
				if(s1 == null) {
					s1 = EMPTY_SET;
				}

				s2 = _mast_results.get(cname);
				if(s2 == null) {
					s2 = EMPTY_SET;
				}

				if((s1.countOnes() + s2.countOnes()) > max_size) {
					best_p1 = p1;
					best_q1 = q1;
					max_size = s1.countOnes() + s2.countOnes();
				}

				// check p1 with q2
				q1 = q1.complement;

				name.set(p1.leaf_vector,q1.leaf_vector);
				cname.set(p1.complement.leaf_vector,q1.complement.leaf_vector);

				s1 = _mast_results.get(name);
				if(s1 == null) {
					s1 = EMPTY_SET;
				}

				s2 = _mast_results.get(cname);
				if(s2 == null) {
					s2 = EMPTY_SET;
				}

				if((s1.countOnes() + s2.countOnes()) > max_size) {
					best_p1 = p1;
					best_q1 = q1;
					max_size = s1.countOnes() + s2.countOnes();
				}
			}
		}

		// build the result
		name.set(best_p1.leaf_vector,best_q1.leaf_vector);
		cname.set(best_p1.complement.leaf_vector,best_q1.complement.leaf_vector);

		BitVector result = new BitVector(_num_leaves);

		s1 = _mast_results.get(name);
		if(s1 != null) {
			result.or(s1);
		}

		s2 = _mast_results.get(cname);
		if(s2 != null) {
			result.or(s2);
		}

		LinkedList<String> leaf_names = new LinkedList<String>();
		int pos = 0;
		for(boolean bit : result) {
			if(bit) {
				leaf_names.add(_lnum2name[pos]);
			}

			pos++;
		}

		tree1.constrainByLeaves(leaf_names);

		return tree1;
	}

	protected void makeLeafSetsEqual(MutableTree t1, MutableTree t2) {

		Iterator<? extends TMutableNode> it;

		it = t1.getNodes().iterator();
		while(it.hasNext()) {
			TMutableNode n = it.next();

			if(n.isLeaf()) {
				if(t2.getNode(n.getName()) == null) {
					n.removeNode();
					it = t1.getNodes().iterator();
				}
			}
		}

		it = t2.getNodes().iterator();
		while(it.hasNext()) {
			TMutableNode n = it.next();

			if(n.isLeaf()) {
				if(t1.getNode(n.getName()) == null) {
					n.removeNode();
					it = t2.getNodes().iterator();
				}
			}
		}

		Trees.removeBinaryNodes(t1);
		Trees.removeBinaryNodes(t2);
	}

	protected boolean leafSetsEquals(Tree t1, Tree t2) {

		for(TNode n : t1.getNodes()) {
			if(n.isLeaf()) {
				if(t2.getNode(n.getName()) == null) {
					return false;
				}
			}
		}

		for(TNode n : t2.getNodes()) {
			if(n.isLeaf()) {
				if(t1.getNode(n.getName()) == null) {
					return false;
				}
			}
		}

		return true;
	}

	private void printSubMASTResult(ESubtree t1, ESubtree t2, BitVector result) {

		// print t1
		System.out.print("{");
		int idx = 0;
		for(Boolean b : t1.leaf_vector) {
			if(b) {
				System.out.print("" + _lnum2name[idx] + " ");
			}

			idx++;
		}
		System.out.print("}; ");

		// print t2
		System.out.print("{");
		idx = 0;
		for(Boolean b : t2.leaf_vector) {
			if(b) {
				System.out.print("" + _lnum2name[idx] + " ");
			}

			idx++;
		}
		System.out.print("} = ");

		// print the result
		System.out.print("{");
		idx = 0;
		for(Boolean b : result) {
			if(b) {
				System.out.print("" + _lnum2name[idx] + " ");
			}

			idx++;
		}
		System.out.print("}");

		System.out.println();
	}

	/**
	 * This method computes the masts for all subtrees in tree1 and tree2 using a dynamic algorithm.  This method
	 * reinitializes all of the internal data structures of SteelWarnowMAST, so it should only be used once per
	 * pair of trees.
	 *
	 * @param tree1 must be unrooted
	 * @param tree2 must be unrooted
	 */
	protected void computeSubMASTs(STITree<ESubtree> tree1, STITree<ESubtree> tree2, boolean is_rooted, boolean cleanup) {

		if(tree1.isRooted() || tree2.isRooted()) {
			throw new RuntimeException("computeSubMASTs: Trees must be unrooted");
		}

		// 0. Assign leaf numbering
		_lname2num = new Hashtable<String,Integer>();
		int next_num = 0;

		for(TNode n : tree1.getNodes()) {
			if(n.isLeaf() && !_lname2num.containsKey(n.getName())) {
				_lname2num.put(n.getName(), next_num++);
			}
		}

		for(TNode n : tree2.getNodes()) {
			if(n.isLeaf() && !_lname2num.containsKey(n.getName())) {
				_lname2num.put(n.getName(), next_num++);
			}
		}

		_num_leaves = next_num;
		_lnum2name = new String[_num_leaves];
		for(Map.Entry<String,Integer> e : _lname2num.entrySet()) {
			_lnum2name[e.getValue()] = e.getKey();
		}

		EMPTY_SET = new BitVector(_num_leaves);

		// 1. find all the e-subtrees for both trees
		_leafset2subtree = new Hashtable<BitVector,ESubtree>();
		_esubtrees1 = generateESubtrees(tree1);
		_esubtrees2 = generateESubtrees(tree2);

		// 2. order them
		orderESubtrees(_esubtrees1);
		orderESubtrees(_esubtrees2);

		// 3. calculate the MAST for each pair
		ESubtreePair rpair = new ESubtreePair();
		_mast_results = new Hashtable<ESubtreePair,BitVector>();

		for(int i = 0; i < _esubtrees1.length; i++) {
			for(int j = 0; j < _esubtrees2.length; j++) {
				BitVector result = calculateMAST(_esubtrees1[i], _esubtrees2[j], _mast_results);

				// only keep informative sets
				if(result.countOnes() > 0) {
					//printSubMASTResult(_esubtrees1[i], _esubtrees2[j], result);

					_mast_results.put(new ESubtreePair(_esubtrees1[i],_esubtrees2[j]), result);
				}

				// delete mast information related to tree2 subtrees that we don't need any more to save space
				if(is_rooted && i == (_esubtrees1.length - 1) && cleanup) {
					for(ESubtree e : _esubtrees2[j].subtrees) {

						if(e.is_root) {
							continue;
						}

						e.decrementRefs();

						if(e.canDeleteRMASTInfo()) {
							for(int k = 0; k < _esubtrees1.length; k++) {
								rpair.set(_esubtrees1[k], e);
								_mast_results.remove(rpair);
							}
						}
					}
				}
			}

			// delete mast information related to tree1 subtrees that we don't need any more to save space.
			if(is_rooted && cleanup) {
				for(ESubtree e : _esubtrees1[i].subtrees) {

					// skip esubtrees that will be included in the final rooted mast computation
					if(e.is_root) {
						continue;
					}

					e.decrementRefs();

					if(e.canDeleteRMASTInfo()) {
						for(int k = 0; k < _esubtrees2.length; k++) {
							rpair.set(e, _esubtrees2[k]);
							_mast_results.remove(rpair);
						}
					}
				}
			}
		}

		return;
	}

	/**
	 * Generates all e-subtrees of the specified tree.
	 *
	 * @param tree is the input tree
	 *
	 * @return all the esubtrees of the input tree
	 */
	protected ESubtree[] generateESubtrees(STITree<ESubtree> tree) {

		ESubtree[] esubtrees = new ESubtree[2 * (tree.getNodeCount() - 1)];

		int next_idx = computeESubtrees(tree.getRoot(), esubtrees, 0);

		assert next_idx == (tree.getNodeCount() - 1);

		computeChildComplements(tree.getRoot(), esubtrees, next_idx);

		return esubtrees;
	}

	/**
	 * This method computes the 'downward' e-subtrees.  This is done by head recursion.
	 *
	 * @return the next index in <code>esubtrees</code> in which an esubtree should be stored.
	 */
	protected int computeESubtrees(STINode<ESubtree> node, ESubtree[] esubtrees, int next_idx) {

		ESubtree[] subtrees = new ESubtree[node.getChildCount()];
		BitVector leaves = new BitVector(_num_leaves);

		if(node.isLeaf()) {
			leaves = getLeafVector(node);
		} else {
			int i = 0;
			for(STINode<ESubtree> child : node.getChildren()) {

				// compute the esubtree for the child
				next_idx = computeESubtrees(child, esubtrees, next_idx);

				// only keep track of child results if we aren't dealing with the root - the root is never involved
				// in these computations.
				if(node != node.getTree().getRoot()) {
					subtrees[i] = child.getData();
					child.getData().incrementRefs();  // increment the number of references to this esubtree
					leaves.or(subtrees[i++].leaf_vector);
				}
			}
		}

		if(node != node.getTree().getRoot()) {
			esubtrees[next_idx] = new ESubtree(this, leaves, subtrees);

			if(node.getParent() == node.getTree().getRoot()) {
				esubtrees[next_idx].is_root = true;
			}

			node.setData(esubtrees[next_idx]);
			return next_idx + 1;
		} else {
			return next_idx;
		}
	}

	/**
	 * This method computes all the complementing esubtrees - those subtrees on the upper part of a branch cut.
	 *
	 * @return the next index in <code>esubtrees</code> in which an esubtree should be stored.
	 */
	protected int computeChildComplements(STINode<ESubtree> node, ESubtree[] esubtrees, int next_idx) {

		// make the complement for each child of this node
		for(STINode<ESubtree> child : node.getChildren()) {
			ESubtree[] subtrees;

			// TODO: Do we need to collapse the root if it has only two children?
			if(node != node.getTree().getRoot()) {
				subtrees = new ESubtree[node.getChildCount()];
			} else {
				// we include one less because the root has no parent
				subtrees = new ESubtree[node.getChildCount() - 1];
			}

			BitVector leaves = new BitVector(_num_leaves);

			int i = 0;
			for(STINode<ESubtree> c : node.getChildren()) {
				if(c != child) {
					subtrees[i++] = c.getData();
					c.getData().incrementRefs();
					leaves.or(c.getData().leaf_vector);
				}
			}

			// add the node-parent complement
			if(node != node.getTree().getRoot()) {
				subtrees[i] = node.getData().complement;
				node.getData().complement.incrementRefs();
				leaves.or(node.getData().complement.leaf_vector);
			}

			ESubtree child_complement = new ESubtree(this, leaves, subtrees);
			child_complement.is_root = child.getData().is_root;
			child.getData().complement = child_complement;
			child_complement.complement = child.getData();
			esubtrees[next_idx++] = child_complement;

			next_idx = computeChildComplements(child, esubtrees, next_idx);
		}

		return next_idx;
	}

	/**
	 * Returns the leaf-vector representation of a subtree.
	 *
	 * @param node is the subtree whose leaf-vector should be constructed
	 *
	 * @return the leaf-vector
	 */
	protected BitVector getLeafVector(TNode node) {

		BitVector leaf_vector = new BitVector(_num_leaves);

		// get leaf names
		setLeafBits(node, leaf_vector);

		return leaf_vector;
	}

	/**
	 * This method is called by {@link getLeafVector} to set the bits in the leaf-vector
	 * that correspond to the leaves underneath the node.
	 *
	 * @param node is the node whose leafs to find
	 * @param leaf_vector is the vector whose bits will be set.
	 */
	protected void setLeafBits(TNode node, BitVector leaf_vector) {

		if(node.isLeaf()) {
			if(!_lname2num.containsKey(node.getName())) {
				//out.printf("Name = %s\n",node.getName());
			}
			leaf_vector.setValue(_lname2num.get(node.getName()), true);
		} else {
			for(TNode child : node.getChildren()) {
				setLeafBits(child, leaf_vector);
			}
		}

		return;
	}

	/**
	 * This method is called by {@link computeSubMASTs} in order to order a list of subtrees by inclusion.
	 * This means that if <code>esubtrees[j]</code> contains the leaves of <code>esubtrees[i]</code>, then
	 * <code>i \lt j</code>.  This method puts the nodes in this order.
	 */
	protected void orderESubtrees(ESubtree[] esubtrees) {

		// reorder the esubtrees according to partial order, then with a linear extension
		// Derek's Note: In this case, the linear extension is equivalent to sorting the
		// e-subtrees by the number of leaves they contain.
		Arrays.sort(esubtrees,new Comparator<ESubtree>() {
			public int compare(ESubtree t1, ESubtree t2) {
				int num_leaves1 = t1.leaf_vector.countOnes();
				int num_leaves2 = t2.leaf_vector.countOnes();

				if(num_leaves1 < num_leaves2) {
					return -1;
				} else if(num_leaves1 == num_leaves2){
					return 0;
				} else {
					return 1;
				}
			}
		});

		// done!
		return;
	}

	/**
	 * This method is called by {@link computeSubMASTs} to find the mast for a pair of e-subtrees.
	 *
	 * @param stree1 is the esubtree from tree 1
	 * @param stree2 is the esubtree from tree 2
	 * @param mast_results a hashtable of the MASTs for pairs of e-subtrees that preceed these
	 * e-subtrees in order (see {@link orderESubtrees}).  The value is a set of the names of leaves in the MAST
	 * of that pair of e-subtrees.
	 *
	 * @return a leaf-vector representation of the solution
	 */
	protected BitVector calculateMAST(ESubtree stree1, ESubtree stree2, Hashtable<ESubtreePair,BitVector> mast_results) {

		// if one of the trees has only one leaf, see if the other subtree contains that leaf
		// if it does, return this one leaf, otherwise return an empty set.
		if(stree1.leaf_vector.countOnes() == 1) {
			BitVector result = new BitVector(stree1.leaf_vector);
			result.and(stree2.leaf_vector);

			return result;
		}

		if(stree2.leaf_vector.countOnes() == 1) {
			BitVector result = new BitVector(stree2.leaf_vector);
			result.and(stree1.leaf_vector);

			return result;
		}

		// otherwise...
		// compute the canonical names for subtrees
		ESubtree[] strees1 = stree1.subtrees;
		ESubtree[] strees2 = stree2.subtrees;

		// calculate the maximal weighted bipartite matching between subtrees of each e-subtree
		BitVector mast = calculateMaxWeightedBipartiteMatch(strees1, strees2, mast_results);
		int mast_size = mast.countOnes();

		// find the max MAST of each e-subtree vs. a subtree of the other e-subtree
		ESubtreePair p = new ESubtreePair();
		for(int i = 0; i < strees2.length; i++) {
			p.set(stree1, strees2[i]);
			BitVector tmp_set = mast_results.get(p);

			if(tmp_set != null && tmp_set.countOnes() > mast_size) {
				mast = tmp_set;
				mast_size = tmp_set.countOnes();
			}
		}

		for(int i = 0; i < strees1.length; i++) {
			p.set(strees1[i], stree2);
			BitVector tmp_set = mast_results.get(p);

			if(tmp_set != null && tmp_set.countOnes() > mast_size) {
				mast = tmp_set;
				mast_size = tmp_set.countOnes();
			}
		}

		return mast;
	}

	/**
	 * This method is called by {@link calculateMAST} to find the maximum match in the weighted bipartite graph
	 * described by the subtrees specified. Currently this method uses an O(n^3) algorithm.
	 *
	 * @param strees1 is the set of leaf-vectors for each of the immediate subtrees in e-subtree 1
	 * @param strees2 is the set of leaf-vectors for each of the immediate subtrees in e-subtree 2
	 * @param mast_results holds the weights for different pairings of subtrees.  The values are leaf-vectors.
	 * Their size is the number of '1's they contain.
	 *
	 * @return the leaf-vector containing the matching.
	 */
	protected BitVector calculateMaxWeightedBipartiteMatch(ESubtree[] strees1, ESubtree[] strees2, Hashtable<ESubtreePair,BitVector> mast_results) {

		// setup the bipartite matcher
		ESubtreePair p = new ESubtreePair();
		HungarianBipartiteMatcher matcher = new HungarianBipartiteMatcher(strees1.length, strees2.length);
		int edges_added = 0;
		for(int i = 0; i < strees1.length; i++) {
			for(int j = 0; j < strees2.length; j++) {
				p.set(strees1[i],strees2[j]);
				BitVector s = mast_results.get(p);

				if(s != null && s.countOnes() > 0) {
					matcher.addEdge(i, j, s.countOnes());
					edges_added++;
				}
			}
		}

		if(edges_added == 0) {
			return EMPTY_SET;
		}

		// find the matching
		matcher.findMatching();

		BitVector result = new BitVector(_num_leaves);

		int[][] matching = matcher.getMatching();

		for(int[] match : matching) {
			p.set(strees1[match[0]],strees2[match[1]]);
			BitVector s = mast_results.get(p);

			result.or(s);
		}

		// Leave this sanity check in.  This catches a surprising number of bugs and has been a huge time-saver
		if(matcher.getMatchingWeight() != Double.NEGATIVE_INFINITY && result.countOnes() != matcher.getMatchingWeight()) {
			// TODO Replace these and other error printouts with proper Exception throwing.
			System.err.println("RESULT " + result.countOnes() + ", doesn't match MatchingWeight " + matcher.getMatchingWeight());
			System.err.println("This may be caused by non-unique leaf-names or leaf-names that contain the '+' character");
		}

		return result;
	}


}
