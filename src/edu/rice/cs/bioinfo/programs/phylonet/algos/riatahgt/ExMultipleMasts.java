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

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/4/11
 * Time: 4:45 PM
 * To change this template use File | Settings | File Templates.
 */

import edu.rice.cs.bioinfo.programs.phylonet.algos.mast.SteelWarnowMAST;
import edu.rice.cs.bioinfo.programs.phylonet.structs.BitVector;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * This class tries to list all possible masts between two trees.
 *
 * @author Cuong Than
 *
 */
public class ExMultipleMasts extends SteelWarnowMAST {
	/* Methods for ExMultipleMasts */
	/**
	 * This method initializes the necessary data members for the computation of multiple MASTs.
	 *     1. First, it calls the base class' function ComputeRMAST so that we can know the MAST for each
	 *        pair of esubtrees. We also get _ordered_esubtrees and _mast_results intialized.
	 *     2. After that, the method fills out the _multipleMastTable. Each element in the table is a pair
	 *        (key, MastStruct), where key = name1 "-" name2, and ExMastStruct stores all information about (e1, e2).
	 *
	 * @return: null if tr1 and tr2 have multiple MASTs; otherwise, return the unique MAST. In case of multiple
	 *     MASTs, the function exComputeMultipleMasts will find all.
	 */
	private Tree setupMultipleMasts(Tree tr1, Tree tr2)
	{
		// After this step, all data members inherited from SteelWarnowMAST are properly initialized.
		Tree mastTree = inner_computeRMAST(tr1, tr2, false);
		_mastTreeSize = mastTree.getLeafCount();
		_multipleMastTable = new Hashtable<ESubtreePair, MastStruct>();

		if (_esubtrees1 == null || _esubtrees2 == null) {
			return mastTree;	// This is the unique MAST.
		}

		// Init the hashtable for all pairs of esubtrees.
		for (int i = 0; i < _esubtrees1.length; i++) {
			ESubtree e1 = _esubtrees1[i];
			for (int j = 0; j < _esubtrees2.length; j++) {
				ESubtree e2 = _esubtrees2[j];
				ESubtreePair key = new ESubtreePair(e1, e2);
				BitVector mast = _mast_results.get(key);

				if (mast != null) {
					int size = mast.countOnes();
					MastStruct ms = new MastStruct(e1, e2, size);
					if (e1.leaf_vector.countOnes() == 1 || e2.leaf_vector.countOnes() == 1) {	// Handle the case of single leaves.
						ms._masts.add(mast);
					}
					_multipleMastTable.put(key, ms);
				}
			}
		}

		_mast_results.clear();	// We don't need to use _mast_results any more.

		return null;			// This pair of tr1 and tr2 might have multiple MASTs.
	}

	/**
	 * Given a list of matchings, this methods enumerates all MASTs.
	 *
	 * @param matchingList
	 *
	 * @return A set of MASTs, where each MAST is a BitVector.
	 */
	private Set<BitVector> getMastSetFromListOfMatching(LinkedList matchingList)
	{
		if (matchingList == null)
			return null;

		Iterator it = matchingList.iterator();
		if (it.hasNext())
			return getMastFromListOfMatchingHelper(it);
		else
			return null;
	}

	/**
	 * The helper function for the above method. This function returns the set of MAST for the subset of the matching list
	 * beginning from the iterator it.
	 *
	 * @param it: Indicates the start of the sublist. When the above function calls with it = matchingList.iterator(),
	 * it obtains the whole set of MASTs.
	 *
	 * @return A set of MASTs for the sublist of the matching list.
	 */
	private Set<BitVector> getMastFromListOfMatchingHelper(Iterator it)
	{
		MastStruct ms = (MastStruct) it.next();
		Set<BitVector> result = null;
		result = new HashSet<BitVector>();

		if (!it.hasNext()) {	// ms is the last element in the list of matchings.
			result = ms._masts;
			return result;
		}

		// Otherwise, there're still more MastStruct in the matching list. Recurse.
		Set<BitVector> subMastSet = getMastFromListOfMatchingHelper(it);

		// and add each element of this Set to each element in the set of sub-Masts of ms.
		for (BitVector leaves : ms._masts) {
			Iterator subMastSetIt = subMastSet.iterator();
			while (subMastSetIt.hasNext()) {
				BitVector temp = null;

				temp = new BitVector(leaves);
				temp.or((BitVector) subMastSetIt.next());	// Grow with new leaves from the remaining pairs.
				result.add(temp);
			}
		}

		return result;
	}

	/**
	 * This method initializes the bipartite graph, which will be used for finding multiple maximum matchings.
	 *
	 * @param names1
	 * @param names2
	 *
	 * @return: The BipartiteGraph
	 */
	private BipartiteGraph setupBipartiteGraph(BitVector names1[], BitVector names2[])
	{
		BipartiteGraph bg = null;
		bg = new BipartiteGraph();
		bg._edgeIndices = new int[names1.length * names2.length][2];
		int edgeCount = 0;

		// Initialize _edgeIndices.
		for(int i = 0; i < names1.length; i++) {
			for(int j = 0; j < names2.length; j++) {
				ESubtreePair key = null;
				key = new ESubtreePair(names1[i], names2[j]);
				MastStruct ms = (MastStruct) _multipleMastTable.get(key);

				if(ms != null && ms._size > 0) {
					bg._edgeIndices[edgeCount][0] = i;
					bg._edgeIndices[edgeCount][1] = j;
					edgeCount++;
				}

				key = null;	// Explicit nulling.
				ms = null;	// Explicit nulling.
			}
		}

		// Setup _edgeBit.
		bg._size = edgeCount;
		bg._edgeBit = new ExBitSet(bg._size);

		return bg;
	}

	/**
	 * This method determines if e1 or e2 of ms is in the list of MastStruct or not.
	 *
	 * @param list
	 * @param ms
	 * @return
	 */
	private boolean isNodeConflict(LinkedList list, MastStruct ms)
	{
		if (list == null)
			return true;
		else {
			Iterator it = list.iterator();
			boolean conflict = false;
			while (it.hasNext() && !conflict) {
				MastStruct temp = (MastStruct) it.next();
				if (ms._e1 == temp._e1 || ms._e2 == temp._e2)
					conflict = true;
			}

			return conflict;
		}
	}

	/**
	 * This function returns a list of pairs of esubtrees (in the form of MastStruct), which forms the maximum bipartite
	 * matching in the graph (stored in bg).
	 *
	 * @param bg
	 * @param names1
	 * @param names2
	 * @param matchingSize
	 *
	 * @return List of MastStruct's if the matching size is really equal to the pre-known matchingSize; null otherwise.
	 */
	private LinkedList<MastStruct> getMatching(BipartiteGraph bg, BitVector names1[], BitVector names2[], int matchingSize)
	{
		LinkedList<MastStruct> matchingList = null;
		matchingList = new LinkedList<MastStruct>();
		int size = 0;

		// Obtain the list of edges in the matching.
		for (int i = 0; i < bg._size; i++) {
			int index1 = bg._edgeIndices[i][0];	// The index in names1.
			int index2 = bg._edgeIndices[i][1];	// The index in names2.

			if (bg._edgeBit.get(i)) {	// This edge is in the matching.
				ESubtreePair key = null;
				key = new ESubtreePair(names1[index1], names2[index2]);
				MastStruct ms = (MastStruct) _multipleMastTable.get(key);
				if (isNodeConflict(matchingList, ms))
					return null;
				else {
					if (ms != null) {
						matchingList.add(ms);
						size += ms._size;
					}
				}

				key = null;	//Explicit nulling.
				ms = null;	//Explicit nulling.
			}
		}

		// We DO obtain a matching. Check with the matchingSize.
		if (size != matchingSize)
			return null;
		else
			return matchingList;
	}

	/**
	 * This method find a set of matchings between two trees whose children's names are in names1 and names2.
	 *
	 * @param names1: The array of immediate children of esubtrees1.
	 * @param names2: The array of immediate children of esubtrees2. esubtree1 and esubtree2 are two trees
	 * that we want to find all maximum matchings.
	 * @param matchingSize: The maximum matching size of esubtree1 and esubtree2 that we obtain from the function
	 * computeRMAST of the base class SteelWarnowMAST.
	 *
	 * @return The set of matchings (in the form of a linked list).
	 */
	private Set<LinkedList<MastStruct>> computeMultipleMatchings(BitVector names1[], BitVector names2[], int matchingSize)
	{
		Set<LinkedList<MastStruct>> multipleMatchingSet = null;
		multipleMatchingSet = new HashSet<LinkedList<MastStruct>>();

		// Initialize the bipartite graph.
		BipartiteGraph bg = setupBipartiteGraph(names1, names2);

		// Find matchings by using the binary counter bg._edgeBit.
		while (bg._edgeBit.increase()) {
			LinkedList<MastStruct> matchingList = getMatching(bg, names1, names2, matchingSize);
			if (matchingList != null)
				multipleMatchingSet.add(matchingList);
		}

		bg._edgeBit.clear();
		bg = null;	// Explicit nulling.
		return multipleMatchingSet;
	}

	/**
	 * This method computes MASTs between the two e-subtrees e1 and e2. All results are updated back to the
	 * data member _multipleMastTable.
	 *
	 * @param e1
	 * @param e2
	 */
	private void computeMultipleMasts(ESubtree e1, ESubtree e2)
	{
		// Get information about the MAST for e1 and e2. (We already know this.)
		ESubtreePair key = null;
		MastStruct ms = null;
		key = new ESubtreePair(e1, e2);
		ms = (MastStruct) _multipleMastTable.get(key);

		if (ms == null || ms._size == 0) {	// MAST is empty; nothing to do.
			return;
		}

		// Get subtrees of e1 and e2.
		BitVector[] names1 = null;
		BitVector[] names2 = null;
		names1 = new BitVector[e1.subtrees.length];
		names2 = new BitVector[e2.subtrees.length];

		for (int i = 0; i < names1.length; i++) {
			names1[i] = e1.subtrees[i].leaf_vector;
		}
		for (int j = 0; j < names2.length; j++) {
			names2[j] = e2.subtrees[j].leaf_vector;
		}

		// First, get the set of MASTs in the form of a matching between e1 and e2.
		Set<LinkedList<MastStruct>> multipleMatchingSet = computeMultipleMatchings(names1, names2, ms._size);
		if (multipleMatchingSet != null && !multipleMatchingSet.isEmpty()) {
			for (LinkedList<MastStruct> matchingList : multipleMatchingSet) {
				Set<BitVector> mastSet = getMastSetFromListOfMatching(matchingList);
				if (mastSet != null && !mastSet.isEmpty()) {
					for (BitVector newMast : mastSet) {
						if (ms._masts == null) {
							ms._masts.add(newMast);
							continue;
						}

						boolean duplicate = false;
						for (BitVector oldMast : ms._masts) {
							BitVector bv = new BitVector(newMast);
							bv.xor(oldMast);
							if (bv.countOnes() == 0) {	// newMast == oldMast.
								duplicate = true;
								break;
							}
						}

						if (!duplicate)
							ms._masts.add(newMast);
					}
				}

				matchingList.clear();
				mastSet = null;			// Explicit nulling.
				matchingList = null;	// Explicit nulling.
			}
		}

		for (LinkedList<MastStruct> ml : multipleMatchingSet) {
			ml.clear();
		}
		multipleMatchingSet.clear();
		multipleMatchingSet = null;		// Explicit nulling.

		// Second, we add the MAST in the form of (e1, subtree-e2) and (subtree-e1, e2) to ms._masts set.
		for (int j = 0; j < names2.length; j++) {
			key = new ESubtreePair(e1.leaf_vector, names2[j]);
			MastStruct temp = (MastStruct) _multipleMastTable.get(key);

			if (temp != null && temp._size == ms._size && !temp._masts.isEmpty()) {
				for (BitVector newMast : temp._masts) {
					if (ms._masts == null) {
						ms._masts.add(newMast);
						continue;
					}
					boolean duplicate = false;
					for (BitVector oldMast : ms._masts) {
						BitVector bv = new BitVector(newMast);
						bv.xor(oldMast);
						if (bv.countOnes() == 0) {	// newMast == oldMast.
							duplicate = true;
							break;
						}
					}

					if (!duplicate) {
						ms._masts.add(newMast);
					}
				}
			}

			key = null;		// Explicit nulling.
			temp = null;	// Explicit nulling.
		}
		for (int i = 0; i < names1.length; i++) {
			key = new ESubtreePair(names1[i], e2.leaf_vector);
			MastStruct temp = (MastStruct) _multipleMastTable.get(key);

			if (temp != null && temp._size == ms._size && !temp._masts.isEmpty()) {
				for (BitVector newMast : temp._masts) {
					if (ms._masts == null) {
						ms._masts.add(newMast);
						continue;
					}
					boolean duplicate = false;
					for (BitVector oldMast : ms._masts) {
						BitVector bv = new BitVector(newMast);
						bv.xor(oldMast);
						if (bv.countOnes() == 0) {	// oldMast == newMast.
							duplicate = true;
							break;
						}
					}

					if (!duplicate) {
						ms._masts.add(newMast);
					}
				}
			}

			key = null;		// Explicit nulling.
			temp = null;	// Explicit nulling.
		}
	}

	/**
	 * This method computes MASTs between the two trees tr1 and tr2.
	 *
	 * @return: A set of MASTs.
	 */
	public Set<Tree> computeMultipleMasts(final Tree tr1, final Tree tr2)
	{
		// Get the base class data members initialized.
		Tree mastTree = setupMultipleMasts(tr1, tr2);
		Set<Tree> result = new HashSet<Tree>();

		if (mastTree != null) {	// We have just one MAST.
			result.add(mastTree);
			return result;
		}

		// We have multiple MASTs here. Fill out the data member _multipleMastTable.
		for (int i = 0; i < _esubtrees1.length; i++) {
			ESubtree e1 = (ESubtree) _esubtrees1[i];

			// Compute multiple sub-MASTs.
			for (int j = 0; j < _esubtrees2.length; j++) {
				ESubtree e2 = (ESubtree) _esubtrees2[j];
				computeMultipleMasts(e1, e2);
			}
		}

		// Find the maximum matchings between xtr1 and xtr2.
		BitVector[] names1 = new BitVector[tr1.getRoot().getChildCount()];
		BitVector[] names2 = new BitVector[tr2.getRoot().getChildCount()];

		Iterator iter = tr1.getRoot().getChildren().iterator();
		for(int i = 0; iter.hasNext(); i++) {
			names1[i] = getLeafVector((TNode) iter.next());
		}
		iter = tr2.getRoot().getChildren().iterator();
		for(int i = 0; iter.hasNext(); i++) {
			names2[i] = getLeafVector((TNode) iter.next());
		}

		Set<LinkedList<MastStruct>> multipleMatchingSet = computeMultipleMatchings(names1, names2, _mastTreeSize);

		for (LinkedList<MastStruct> matchingList : multipleMatchingSet) {
			Set<BitVector> temp = getMastSetFromListOfMatching(matchingList);	// Set of MASTs.
			if (temp == null) {
				continue;
			}
			for (BitVector mast : temp) {
				mastTree = getTreeFromBits(tr1, mast);
				if (mastTree != null && mastTree.getRoot().getID() == tr1.getRoot().getID())
					result.add(mastTree);
			}
		}

		_multipleMastTable.clear();
		return result;
	}


	/**
	 * This method implements the command-line interface for this class.
	 */


	/**
	 * This method gets the tree whose leaf set is represented by a bit vector. The tree is constrained
	 * on tr.
	 *
	 * @param tr: The tree on which the new tree is constrained.
	 * @param bv: The bit vector reprensenting which leaf is to be maintained.
	 *
	 * @return: The STITree corresponding to bv.
	 */
	private Tree getTreeFromBits(Tree tr, BitVector bv)
	{
		LinkedList<String> leafNames = new LinkedList<String>();
		int pos = 0;

		for (boolean bit : bv) {
			if (bit) {
				leafNames.add(_lnum2name[pos]);
			}
			pos++;
		}

		STITree<Object> result = new STITree<Object>(tr.getRoot(), true);
		result.constrainByLeaves(leafNames);

		return result;
	}

	/**
	 * This method collapses the two trees tr1 and tr2. In the present version, we only implement the
	 * first rule in a paper by Steel. For rule 2, we are still unclear about how to deal with non-binary trees.
	 *
	 * @param tr1: The
	 * @param tr2
	 */
	public void collapse(Tree tr1, Tree tr2)
	{
		runRule1(tr1.getRoot(), tr2);
	}

	/**
	 * The method collapses a subtree, rooted at peerNode1, and tr2. To refine the tree tr1 and tr2,
	 * calls this function with peerNode1 as the root of tr1.
	 *
	 * @param peerNode1
	 * @param tr2
	 *
	 * @return: true if there's actually a collapse; false if there's not.
	 */
	private boolean runRule1(TNode peerNode1, Tree tr2)
	{
		// Base case: nothing to do with a single leaf.
		if (peerNode1.isLeaf())
			return false;

		// Non-trivial case. Recurse on peerNode1's childrent before collapse itself.
		boolean changed = false;

		// Collapse peerNode1's children.
		Iterator it = peerNode1.getChildren().iterator();
		TMutableNode currentNode = null;

		while (it.hasNext()) {
			TMutableNode oldNode = currentNode;
			currentNode =  (TMutableNode) it.next();

			if (oldNode != null) {
				boolean tmp = runRule1(oldNode, tr2);
				changed = changed || tmp;
			}
		}
		if (currentNode != null) {
			boolean tmp = runRule1(currentNode, tr2);
			changed = changed || tmp;
		}

		// Then, test if we can make the whole tr1 and tr2 as single leaves.
		STITree subtree1 = new STITree(peerNode1, true);
		if (isPendant(subtree1)) {
			TMutableNode peerNode2 = findPeerNode(tr2, subtree1);
			if (peerNode2 != null) {
				String newName = peerNode1.getName();
				if (newName == TNode.NO_NAME)
					createNewNameFromLeaves(peerNode2);

				fixNode((TMutableNode) peerNode1, newName);	// Collapse in tr1.
				fixNode((TMutableNode) peerNode2, newName);	// Collapse in tr2.

				changed = true;
			}
		}

		return changed;
	}

	/**
	 * Tests if this tree is pendant, i.e. its immediate children are leaves.
	 *
	 * Require:
	 *     Tree tr is rooted.
	 *
	 * @param tr: The tree to be checked if it's pendant.
	 *
	 * @return
	 *     true if this tree is pedant; false if it's not.
	 */
	private boolean isPendant(Tree tr)
	{
		Iterator it = tr.getRoot().getChildren().iterator();
		boolean pendant = true;

		// Test if all of the tree's children are leaves.
		while (it.hasNext() && pendant) {
			TNode node = (TNode) it.next();

			if (!node.isLeaf())
				pendant = false;
		}

		return pendant;
	}

	/**
	 * Finds a peer node in the tree tr whose subtree is isormorphic to pendantStr.
	 * Require:
	 *     As this function is used by the collapse function, we only provide the implemetation
	 *     in the case where pendantStr is a pendant subtree.
	 * @param
	 *     pendantStr: A pendant subtree, i.e. its immediate children are leaves.
	 * @return
	 *     a node in this tree if pendantStr is actually its subtree; null otherwise.
	 *
	 * NOTE: This method is used in collpasing two trees.
	 */
	private TMutableNode findPeerNode(Tree tr, Tree pendantSubtree)
	{
		// Get the name of first child of pendantStr.
		Iterator it = pendantSubtree.getRoot().getChildren().iterator();
		String name = ((TNode) it.next()).getName();

		// Get the parent of this node in this tree.
		TMutableNode parent = (TMutableNode) tr.getNode(name).getParent();

		// Compare the pendant tree with the subtree rooted at node parent.
		STITree peerSubtree = new STITree(parent, true);
		if (isIdentical(pendantSubtree, peerSubtree))
			return parent;
		else
			return null;
	}

	/**
	 * This method checks if the two trees tr1 and tr2 are identical. For the purpose of
	 * collapsing two trees, this function will check ONLY pendant trees. It'll always
	 * return false if either tr1 and tr2 is not pedant, even if they are really isomorphic.
	 *
	 * @param
	 *     tr1: Pendant tree 1.
	 * @param
	 *     tr2: Pendant tree 2.
	 * @return
	 *     true if tr1 and tr2 are pendant and isomorphic; false otherwise.
	 */
	private boolean isIdentical(Tree tr1, Tree tr2)
	{
		// Make sure that we only compare two trees with the same number of leaves and
		// the same number of internal nodes.
		if (tr1.getLeafCount() != tr2.getLeafCount() || tr1.getNodeCount() != tr2.getNodeCount())
			return false;

		// Forms the array of names of the root's children.
		TNode root = tr1.getRoot();
		String childNames[] = new String[root.getChildCount()];
		Iterator it = root.getChildren().iterator();

		for (int i = 0; it.hasNext(); i++) {
			TNode child = (TNode) it.next();
			if (!child.isLeaf())	// tr1 is not pendant.
				return false;
			else
				childNames[i] = child.getName();
		}

		// Is every child of pendantStr in this tree?
		boolean identical = true;
		it = tr2.getRoot().getChildren().iterator();

		while (it.hasNext() && identical) {
			TNode child = (TNode) it.next();
			String name = child.getName();
			int j;

			if (!child.isLeaf()) {	// tr2 is pendant.
				return false;
			}
			else {
				// Find this child in tr1's children.
				for (j = 0; j < childNames.length; j++) {
					if (name.equals(childNames[j]))
						break;
				}
				if (j >= childNames.length)
					identical = false;
			}
		}

		return identical;
	}

	/**
	 * Creates the new name for the TNode node. The new name is the concatenation of (sorted)
	 * names of node's children, or its current name if it is not an empty string.
	 * Require:
	 *     The subtree rooted at node must be a pendant tree.
	 * @return
	 *     The new name of the node.
	 */
	private String createNewNameFromLeaves(TNode node)
	{
		// Get the names of its children.
		String childNames[] = new String[node.getChildCount()];
		Iterator it = node.getChildren().iterator();

		for (int i = 0; it.hasNext(); i++) {
			TNode child = (TNode) it.next();
			childNames[i] = child.getName();
		}

		// Concatenate names from childNames to create a new name for the TNode node.
		String newName = new String();
		Arrays.sort(childNames);
		for (int i = 0; i < childNames.length; i++) {
			if (i > 0) {
				newName += "-";
			}
			newName += childNames[i];
		}

		return newName;
	}

	/**
	 * This method actually deletes all children of the TNode node, and renames its name
	 * as the concatenation of the names of the children. Note that the names are sorted
	 * alphabetically before they are concatenated.
	 *
	 * Require:
	 *     Again, node must be the the root of a pedant subtree.
	 * @param
	 *     node: The nod to be fixed.
	 * @param
	 *     newName: Name for the new node.
	 *
	 * NOTE: This method is used in the collapse procedure.
	 */
	private void fixNode(TMutableNode node, String newName)
	{
		// Remove children of this node.
		Iterator it = node.getChildren().iterator();
		while (it.hasNext()) {
			TMutableNode child = (TMutableNode) it.next();
			node.removeChild(child, false);
			it = node.getChildren().iterator();
		}

		// Set the node's name as newName.
		String oldName = node.getName();
		if (oldName != newName)
			node.setName(newName);
	}

	/**
	 * This method refines two trees tr1 and tr2.
	 *
	 * @param tr1: (Mutable) tree 1.
	 * @param tr2: (Mutable) tree 2.
	 *
	 * @return
	 *     true if there's actually refinment on tr1 and tr2; false otherwise.
	 */
	public void refine(Tree tr1, Tree tr2)
	{
		// Generate all possible splits.
		LinkedList<Set<String>> splitList1 = new LinkedList<Set<String>>();
		LinkedList<Set<String>> splitList2 = new LinkedList<Set<String>>();

		generateAllSplits(tr1, splitList1);
		generateAllSplits(tr2, splitList2);

		// Find compatible splits on tr1 and refine tr2 on them.
		Set<String> taxa = getLeafSet(tr1);
		Iterator<Set<String>> it;

		it = splitList1.iterator();
		while (it.hasNext()) {
			Set<String> split1 = it.next();

			if (isCompatible(split1, splitList2, taxa)) {
				Set<String> temp = new HashSet<String>(split1);
				TMutableNode refinementNode = (TMutableNode) findPeerNode(tr2, temp, taxa);
				fixNode(tr2, refinementNode, temp);
			}
		}

		// Find compatible splits on tr2 and refine tr1 on them.
		it = splitList2.iterator();
		while (it.hasNext()) {
			Set<String> split2 = (Set<String>) it.next();

			if (isCompatible(split2, splitList1, taxa)) {
				Set<String> temp = new HashSet<String>(split2);
				TMutableNode refinementNode = (TMutableNode) findPeerNode(tr1, temp, taxa);
				fixNode(tr1, refinementNode, temp);
			}
		}
	}

	/**
	 * The function generates all splits for a tree.
	 *
	 * @param
	 *     tr: The tree that the function returns its splits.
	 * @param
	 *     splitList: Holds the splits.
	 * Require:
	 *     splitList must be allocated memory before calling this function.
	 */
	private void generateAllSplits(Tree tr, LinkedList<Set<String>> splitList)
	{
		// Base case. There're no splits for a tree with a single leaf.
		if (tr.getNodeCount() == 1)
			return;

		// General case.
		Iterator it = tr.getRoot().getChildren().iterator();
		while (it.hasNext()) {
			STITree subtree = new STITree((TNode) it.next(), true);
			Set<String> leafSet = getLeafSet(subtree);

			if (leafSet.size() > 1) {
				splitList.add(leafSet);
				generateAllSplits(subtree, splitList);
			}
		}
	}


	/**
	 * The method returns the set of leaves of the tree tr.
	 *
	 * @param
	 *     tr: The input tree
	 * @return
	 *     The set of names of leaves.
	 */
	private Set<String> getLeafSet(Tree tr)
	{
		Set<String> leafSet = new HashSet<String>();

		// Single leaf.
		if (tr.getRoot().isLeaf()) {
			leafSet.add(tr.getRoot().getName());
			return leafSet;
		}

		// General tree.
		Iterator it = tr.getRoot().getChildren().iterator();
		while (it.hasNext()) {
			TNode node = (TNode) it.next();
			leafSet.addAll(getLeafSet(node));
		}

		return leafSet;
	}

	/**
	 * The method performs the same functionality as the function getLeafSet(tree), but
	 * it now accepts a node as a parameter.
	 *
	 * @param
	 *     node: The root of the subtree that the funtion gets its leaves.
	 * @return
	 *     Set of leaves of the subtree rooted at the TNode node.
	 */
	private Set<String> getLeafSet(TNode node)
	{
		Tree tr = new STITree(node, true);
		return getLeafSet(tr);
	}

	/**
	 * This method tests if the bipartition split1 is compatible with all bipartitions in the splitList.
	 *
	 * @param
	 *     split1: The split to be checked for compatibility.
	 * @param
	 *     splitList: The list of splits to be checked against split1.
	 * @param
	 *     taxa: The set of names of taxa.
	 * @return
	 *     true if split1 is compatible with, all splits in splitList; false otherwise.
	 */
	private boolean isCompatible(Set<String> split1, LinkedList<Set<String>> splitList, Set<String> taxa)
	{
		Iterator<Set<String>> it = splitList.iterator();
		boolean compatible = true;

		// Test split1 for compatibility with every split in splitList.
		while (it.hasNext() && compatible) {
			Set<String> split2 = it.next();
			if (!isCompatibleHelper(split1, split2, taxa)) {
				return false;
			}
		}

		return compatible;
	}

	/**
	 * This function tests if two bipartitions A1|B1 and A2|B2 are compatible.
	 *
	 * @param
	 *     split1: The set of names of leaves in split 1.
	 * @param
	 *     split2: The set of names of leaves in split 2.
	 * @param
	 *     taxa: The set of names of taxa.
	 * @return
	 *     true if split1 and split2 is compatible; false otherwise.
	 * Require:
	 *     split1 and split2 must not be identical. If they are identical, they are certainly
	 *     compatible and the check for their equality is done before calling this function.
	 */
	private boolean isCompatibleHelper(Set<String> split1, Set<String> split2, Set<String> taxa)
	{
		// Get the two complements of split1 and split2.
		Set<String> complement1 = new HashSet<String>(taxa);
		Set<String> complement2 = new HashSet<String>(taxa);
		complement1.removeAll(split1);
		complement2.removeAll(split2);

		// Check for compatibility.
		Set<String> temp = new HashSet<String>();

		temp.addAll(split1);
		temp.retainAll(split2);
		if (temp.isEmpty())		// A1 \cup A2 = \emptyset.
			return true;

		temp.clear();
		temp.addAll(split1);
		temp.retainAll(complement2);
		if (temp.isEmpty())		// A1 \cup B2 = \emptyset.
			return true;

		temp.clear();
		temp.addAll(complement1);
		temp.retainAll(split2);
		if (temp.isEmpty())		// B1 \cup A2 = \emptyset.
			return true;

		temp.clear();
		temp.addAll(complement1);
		temp.retainAll(complement2);
		if (temp.isEmpty())		// B1 \cup B2 = \emptyset.
			return true;

		// All of the four intersections are not empty. The two splits are not compatible.
		return false;
	}

	/**
	 * This method finds the (unique) node v in the tree tr where we perform the refinement.
	 *
	 * @param
	 *     tr: The tree on which we'll perform the refinement.
	 * @param
	 *     split: A compatible split.
	 * @param
	 *     taxa: The set of names of leaves in the tree tr.
	 * @return
	 *     The node v that corresponds to the compatible split.
	 *
	 * NOTE: This method is used for refinement.
	 */
	private TNode findPeerNode(Tree tr, Set<String> split, Set<String> taxa)
	{
		TNode lca1 = getLca(tr, split);
		if (!lca1.isRoot()) {	// LCA(split) != root.
			Set<String> ls1 = getLeafSet(lca1);
			if (ls1.equals(split))
				return null;	// No refinement.
			else
				return lca1;	// There's a refinement.
		}
		else {	// LCA(split) == root.
			Set<String>	complement = new HashSet<String>(taxa);
			complement.removeAll(split);
			TNode lca2 = getLca(tr, complement);

			if (!lca2.isRoot()) {
				Set<String> ls2 = getLeafSet(lca2);
				if (!ls2.equals(complement)) {
					// We actually refine at lca2.
					split.clear();
					split.addAll(complement);	// Correct split under the refinement node.
					return lca2;
				}
				else {
					return lca2.getParent();
				}
			}
			else {	// Both lca1 and lca2 are the root.
				return lca1;
			}
		}
	}

	private TNode getLca(Tree tr, Set<String> leafSet)
	{
		Iterator<String> it = leafSet.iterator();
		String name = (String) it.next();
		ArrayList<TNode> path1 = getPath(tr.getNode(name));
		TNode lca = tr.getNode(name);

		while (it.hasNext()) {
			name = it.next();
			ArrayList<TNode> path2 = getPath(tr.getNode(name));

			lca = getLca(path1, path2, lca);
		}

		return lca;
	}

	/**
	 * This method returns the path from the TNode node upto the root.
	 *
	 * @param
	 *     node: The starting node.
	 * @return
	 *     The array of nodes on the path.
	 */
	private ArrayList<TNode> getPath(TNode node)
	{
		ArrayList<TNode> path = new ArrayList<TNode>();
		TNode temp = node;

		do {
			path.add(temp);
			temp = temp.getParent();
		} while (temp != null);

		return path;
	}

	/**
	 * Return the least common node in the two paths.
	 * Require:
	 *     path1 and path2 always share a node, at least the root. We start there and goes
	 *     down to find the least commond node.
	 * @param
	 *     path1:
	 * @param
	 *     path2
	 * @param lastLca: the node indicating the node that lca(path1, path2) must be
	 *     its ancestors. When computing lca for multiple leaves, lastLca is the lca of
	 *     all previous nodes.
	 * @return
	 *     The least common ancestor.
	 */
	private TNode getLca(ArrayList<TNode> path1, ArrayList<TNode> path2, TNode lastLca)
	{
		// Find the least common ancestor.
		int i = path1.size() - 1;
		int j = path2.size() - 1;
		TNode lca = null;
		boolean stop = false;
		boolean differ = false; // See if two paths begin to differ.

		while ((i >= 0 && j >= 0) && !stop) {
			lca = (TNode) path1.get(i);
			if (lca.equals(lastLca)) {
				stop = true;
			}
			if (!lca.equals((TNode) path2.get(j))) {
				stop = true;
				differ = true;
			}
			else {
				i--;
				j--;
			}
		}

		if (lca.isRoot())
			return lca;
		else if (differ)
			return lca.getParent();
		else
			return lca;
	}

	/**
	 * This function fixes the node v.
	 *
	 * @param
	 *     tr: The tree to be refined.
	 * @param
	 *     v: The node to be fixed. It's returned by findPeerNode(Tree, split, taxa).
	 * @param
	 *     split: The split that this function uses to fix node v.
	 * @return
	 *     false if we do not modify the tree; true otherwise.
	 * Require:
	 *     split must always be a subset of the the set of leaves under the node v. This
	 *     requirement is done when calling findNode.
	 *
	 * NOTE: This function is used to refine two trees.
	 */
	private boolean fixNode(Tree tr, TMutableNode v, Set<String> split)
	{
		if (v == null) {
			return false;	// No refinement.
		}

		// Fix the node v by breaking it into two distinct u1 (= v) and u2, so that every subtree
		// in split is attached to u1 while the other ones are added under u2.
		STINode parent = (STINode) v.getParent();
		TMutableNode u2;
		if (parent != null)	{	// v is an internal node.
			u2 = parent.createChild();
			Iterator it = v.getChildren().iterator();

			while (it.hasNext()) {
				TMutableNode child = (TMutableNode) it.next();
				Set leafSubset = this.getLeafSet(child);
				leafSubset.retainAll(split);
				if (leafSubset.isEmpty()) {	// This subset of leaves is not in split, so move it under u2.
					u2.adoptChild(child);
					it = v.getChildren().iterator();
				}
			}
			u2.adoptChild(v);

			// Check the special case when v has only one child.
			if (v.getChildCount() == 1) {
				it = v.getChildren().iterator();
				u2.adoptChild((TMutableNode) (it.next()));	// Adopt the only child of v to u2.
				v.removeNode();
			}

			return true;
		}
		else {	// v is the root.
			u2 = v.createChild();

			// Perform the same refinement as above, except that we now move subtrees in split under u2 instead of v as above.
			Iterator it = v.getChildren().iterator();
			while (it.hasNext()) {
				TMutableNode child = (TMutableNode) it.next();

				if (child.getID() == u2.getID()) {
					continue;
				}
				Set leafSubset = getLeafSet(child);
				if (split.containsAll(leafSubset)) {
					u2.adoptChild(child);
					it = v.getChildren().iterator();
				}
			}

			// Check the special case when u2 has only one child.
			if (u2.getChildCount() == 1) {
				it = u2.getChildren().iterator();
				v.adoptChild((TMutableNode) (it.next()));
				u2.removeNode();
			}

			return true;
		}
	}

	/* Data members for ExMultipleMasts */
	private Hashtable<ESubtreePair, MastStruct> _multipleMastTable;	// The table (canonical name, ExMastStruct) for all pairs of esubtrees.
	private int _mastTreeSize;	// The size of the final MAST for the two tree tr1 and tr2.

	/* Classes used internally by ExMultipleMasts */
	class MastStruct {
		/* Methods for MastStruct */
		MastStruct(ESubtree e1, ESubtree e2, int size)
		{
			_e1 = e1;
			_e2 = e2;
			_size = size;
			_masts = new HashSet<BitVector>();
		}

		MastStruct()
		{
			_e1 = null;
			_e2 = null;
			_size = 0;
			_masts = new HashSet<BitVector>();
		}

		/* Data members for MastStruct */
		public ESubtree _e1;			// e-subtree 1.
		public ESubtree _e2;			// e-subtree 2. NOTE: These both data members are for debug only.
		public Set<BitVector> _masts;	// The Set of MASTs. Each MAST is represented as a BitVector.
		public int _size;				// Size of the MAST.
	}

	class ExBitSet extends BitSet {
		ExBitSet(int size)
		{
			super(size);
			clear();
			_bitSetLength = size;
		}

		/**
		 * This method increases the bitset by 1.
		 *
		 * @return: true if we can still increment the bit set, that is, we haven't reached 11...1;
		 *     false if we increment bit set 11...1.
		 */
		public boolean increase()
		{
			int i;

			for (i = 0; i < _bitSetLength && get(i); i++) {
				clear(i);
			}
			if (i < _bitSetLength) {
				set(i);
				return true;	// We can still increase.
			}
			else
				return false;	// We reach the upper limit.
		}

		// Data members.
		int _bitSetLength;
		static final long serialVersionUID = 1;
	}

	class BipartiteGraph {
		public int _edgeIndices[][];	// [i][0] = index of the node in names1. [i][1] = index of node in names2.
		public ExBitSet _edgeBit;		// For finding subsets.
		public int _size;				// The number of edges in the graph.
	}
}

