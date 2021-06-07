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

package edu.rice.cs.bioinfo.programs.phylonet.algos.fitchpars;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.Arrays;
import java.util.Hashtable;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/31/11
 * Time: 3:35 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParsimonyCalculator
{
    // constants
	protected static int MISSING_INT = -2;

	// inner classes
	/**
	 * This class holds parsimony-computation related information for each node in the tree.
	 */
	protected static class ParsimonyInfo {

		// fields
		/**
		 * The sequence at this node.  Generated for internal nodes, set in leaves.
		 */
		public int[] sequence;

		/**
		 * A data structure to hold the list of possible values that can be assigned to this node, when
		 * considering the parsimony score of its subtree.
		 */
		public LinkedList<Integer> assignment = new LinkedList<Integer>();

		/**
		 * The node in the input tree that this node corresponds to.
		 */
		public TMutableNode node;

		// constructors
		public ParsimonyInfo(TMutableNode node) {
			this.node = node;
		}
	}

	// fields
	// this field is used to convert int sequences back into strings
	protected char[] _int2char = new char[256];

	protected char _missing_symbol = ' ';

	// methods
	/**
	 * Set the character that will be treated as a missing or gap character.  By default, the gap character is a space ' '
	 * character.
	 */
	public void setMissingCharacter(char c) {
		_missing_symbol = c;
	}

	/**
	 * @return the character that will be interpretted to be a gap or missing character in the data.
	 */
	public char getMissingCharacter() {
		return _missing_symbol;
	}

	/**
	 * Compute the parsimony score of the tree given the specified set of
	 * sequences.  The branch lengths will be set to be the number of mutations
	 * along that edge.
	 *
	 * @param tree is the tree to compute the parsimony score for
	 * @param taxa_names is the set of leaf names in the tree
	 * @param sequences is the set of sequences.  The sequences will
	 * be mapped to leaves of the tree by mapping the index of each sequence
	 * to the index of a taxa name in the <code>taxa_names</code> parameters
	 *
	 * @return the parsimony score of the tree.
	 */
	public int computeParsimony(MutableTree tree, String[] taxa_names, String[] sequences) {

		return computeParsimony(tree, taxa_names, sequences, null);
	}

	/**
	 * This method computes the parsimony score of the tree and also returns the assignment of sequences to the internal nodes.
	 *
	 * @param assignments is a hashtable that will be populated with the mapping between nodes in the input tree and
	 * the sequences that were assigned to them during this computation.
	 */
	public int computeParsimony(MutableTree tree, String[] taxa_names, String[] sequences, Hashtable<TNode,String> assignments) {

		// convert sequences to integers
		int[][] iseqs = new int[sequences.length][sequences[0].length()];;
		int max_val = convertStrSeqs2Ints(sequences, iseqs);

		// make the hashtable lookup
		assert taxa_names.length == sequences.length;

		Hashtable<String,int[]> names2seqs = new Hashtable<String,int[]>();
		for(int i = 0; i < taxa_names.length; i++) {
			names2seqs.put(taxa_names[i], iseqs[i]);
		}

		// make the tree that will hold the assignments
		STITree<ParsimonyInfo> ptree = new STITree<ParsimonyInfo>();

		// run parsimony computation
		int pscore = computeParsimony(tree, names2seqs, max_val, sequences[0].length(), ptree);

		// record the assignments
		if(assignments != null) {
			for(STINode<ParsimonyInfo> node : ptree.getNodes()) {
				assignments.put(node.getData().node, convertIntSeq2Str(node.getData().sequence));
			}
		}

		return pscore;
	}

	/**
	 * Convert an integer sequence back into a character string.
	 * @param seq
	 * @return
	 */
	protected String convertIntSeq2Str(int[] seq) {

		char[] cseq = new char[seq.length];

		for(int i = 0; i < seq.length; i++) {
			if(seq[i] == MISSING_INT) {
				cseq[i] = _missing_symbol;
			} else {
				cseq[i] = _int2char[seq[i]];
			}
		}

		return new String(cseq);
	}

	/**
	 * Convert a set of string sequences into a set of integer sequences
	 * such that each distinct character symbol receives a distinct integer.
	 *
	 * @param sequences is the set of string sequences
	 *
	 * @return is the integer sequences that represent the string sequences.
	 */
	protected int convertStrSeqs2Ints(String[] sequences, int[][] seqs) {

		int[] assignments = new int[256];
		int next_assignment = 0;

		Arrays.fill(assignments, -1);

		// set the blank assignment
		assignments[_missing_symbol] = MISSING_INT;

		// initialize the int seqs
		for(int i = 0; i < seqs.length; i++) {

			// sanity check
			assert sequences[i].length() == sequences[0].length();

			for(int j = 0; j < seqs[i].length; j++) {
				char c = sequences[i].charAt(j);
				int c_idx = (int) c;

				if(assignments[c_idx] == -1) {
					_int2char[next_assignment] = c;
					assignments[c_idx] = next_assignment++;
				}

				seqs[i][j] = assignments[c_idx];
			}
		}

		return next_assignment;
	}

	/**
	 * Do the actual computation of the parsimony score for the tree and sequences.
	 *
	 * @param max_val is the maximum integer value that appears in the sequence arrays.
	 */
	protected int computeParsimony(MutableTree tree, Hashtable<String,int[]> names2seqs, int max_val, int seq_len, STITree<ParsimonyInfo> ptree) {

		// copy the tree
		copyNode(ptree.getRoot(), tree.getRoot(), names2seqs, seq_len);

		for(int pos = 0; pos < seq_len; pos++) {
			// make initial assignments
			makeInitialAssignments(pos, ptree.getRoot(), max_val);

			// make final value selection
			makeFinalSelection(pos, ptree.getRoot(), ptree.getRoot().getData().assignment.getFirst());
		}

		// compute scores
		return computePScore(ptree.getRoot());
	}

	/**
	 * Compute the parsimony score of the subtree pnode.  This method also updates the branch lengths of this subtree (and the node's incoming branch)
	 * to be the number of mutations along that branch.  This method *also* updates the branch lengths in the original
	 * tree provided as input.
	 */
	protected int computePScore(STINode<ParsimonyInfo> pnode) {

		int pscore = 0;

		// set the parent distance to be the number of mutations between this node and its parent
		if(pnode.isRoot()) {
			pnode.setParentDistance(0);
		} else {
			int[] ps = pnode.getParent().getData().sequence;
			int[] seq = pnode.getData().sequence;

			int num_diffs = 0;

			for(int i = 0; i < ps.length; i++) {
				if(ps[i] != seq[i]) {
					num_diffs++;
				}
			}

			pnode.setParentDistance(num_diffs);
		}

		// set the mirror node's distance
		pnode.getData().node.setParentDistance(pnode.getParentDistance());

		// compute all children's pscores
		for(STINode<ParsimonyInfo> child : pnode.getChildren()) {
			pscore += computePScore(child) + child.getParentDistance();
		}

		return pscore;
	}

	/**
	 * Compute the possible assignments of characters to node <code>pnode</code> considering its subtree for
	 * position <code>pos</code> in the input positions.
	 *
	 * @param max_val is the maximum value that the integer (character) will have.
	 */
	protected void makeInitialAssignments(int pos, STINode<ParsimonyInfo> pnode, int max_val) {
		ParsimonyInfo pi = pnode.getData();

		if(pnode.isLeaf()) {
			pi.assignment.clear();

			// by not making any record of the missing character, we don't bias the choice of ancestral
			// state at all.
			if(pi.sequence[pos] != MISSING_INT) {
				pi.assignment.addFirst(pi.sequence[pos]);
			}
		} else {
			int[] count = new int[max_val];

			// compute all children assignments
			for(STINode<ParsimonyInfo> child : pnode.getChildren()) {
				makeInitialAssignments(pos, child, max_val);

				for(Integer a : child.getData().assignment) {
					count[a]++;
				}
			}

			// find the max and set the initial assignment
			int max = 0;
			pi.assignment.clear();
			for(int i = 0; i < count.length; i++) {
				if(count[i] > max) {
					max = count[i];
					pi.assignment.clear();
					pi.assignment.add(i);
				} else if(count[i] == max) {
					pi.assignment.add(i);
				}
			}
		}

		return;
	}

	/**
	 * After having found the set of possible characters for the subtree, <code>pnode</code>, this
	 * method computes the character that yields the parsimony score of the tree for position <code>pos</code>
	 * in the input sequences.
	 *
	 * @param passignment is the final assignment made to the parent of this node.
	 */
	protected void makeFinalSelection(int pos, STINode<ParsimonyInfo> pnode, int passignment) {
		ParsimonyInfo pi = pnode.getData();
		int a = -1;

		if(pi.assignment.size() > 1) {

			if(pi.assignment.contains(passignment)) {
				a = passignment;
			} else {
				a = pi.assignment.getFirst();
			}
		} else if(pi.assignment.size() == 1) {
			a = pi.assignment.getFirst();
		} else { // the node is entirely unspecified, just use it's parent's node
			a = passignment;
		}

		// set the assignment
		pi.sequence[pos] = a;

		// propagate down to children
		for(STINode<ParsimonyInfo> child : pnode.getChildren()) {
			makeFinalSelection(pos, child, a);
		}

		return;
	}

	/**
	 * Copy the node and all its children.  This method is for creating a tree that is topologically identical
	 * to the input tree.  The resulting tree has additional fields for holding data used to compute the
	 * parsimony score of the tree.
	 */
	protected void copyNode(STINode<ParsimonyInfo> pnode, TMutableNode n, Hashtable<String,int[]> names2seqs, int seq_len) {
		ParsimonyInfo pi = new ParsimonyInfo(n);

		pnode.setData(pi);

		if(n.isLeaf()) {
			pi.sequence = names2seqs.get(n.getName());

			if(pi.sequence == null) {
				throw new RuntimeException("No sequence provided for leaf node " + n.getName());
			}
		} else {
			pi.sequence = new int[seq_len];
			for(TMutableNode c : n.getChildren()) {
				STINode<ParsimonyInfo> pchild = pnode.createChild();
				copyNode(pchild, c, names2seqs, seq_len);
			}
		}

		return;
	}
}
