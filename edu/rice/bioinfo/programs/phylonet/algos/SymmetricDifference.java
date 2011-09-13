package edu.rice.bioinfo.programs.phylonet.algos;

import edu.rice.bioinfo.programs.phylonet.structs.BitVector;
import edu.rice.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.bioinfo.programs.phylonet.structs.tree.util.Bipartitions;
import edu.rice.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Hashtable;
import java.util.Map;

/**
 * This is an implementation of the Symmetric Difference tree metric.  It computes the number of
 * bipartitions in tree 2 that are not in tree 1 (false negatives), and the number of bipartitions
 * in tree 1 that are not in tree 2 (false positives).
 *
 * @author Derek Ruths
 */
public class SymmetricDifference {

	// fields
	protected int _fp;
	protected int _fn;

	protected double _wfp;
	protected double _wfn;

	protected int _num_iedges1 = -1;
	protected int _num_iedges2 = -1;

	// methods
	/**
	 * Compute the symmetric difference between two trees.  The false positives
	 * and false negatives can be retrieved through the methods
	 * {@link getFalsePositives} and {@link getFalseNegatives}, respectively.
	 *
	 * @throws RuntimeException if the trees do not have identical leaf sets
	 */
	public void computeDifference(Tree t1, Tree t2) {

		if(!Trees.leafSetsAgree(t1, t2)) {
			throw new RuntimeException("Trees must have identical leaf sets");
		}

		// compute bipartitions
		Map<String,Integer> la = Bipartitions.assignLeafPositions(t1);
		Map<BitVector,TNode> bipartitions1 = new Hashtable<BitVector,TNode>();
		Bipartitions.computeBipartitions(t1, la, bipartitions1);

		Map<BitVector,TNode> bipartitions2 = new Hashtable<BitVector,TNode>();
		Bipartitions.computeBipartitions(t2, la, bipartitions2);

		// count missing bipartitions
		_fp = 0;
		for(BitVector bv : bipartitions1.keySet()) {
			if(!bipartitions2.containsKey(bv)) {
				bv.not();
				if(!bipartitions2.containsKey(bv)) {
					_fp++;
				}
				bv.not();
			}
		}

		_fn = 0;
		for(BitVector bv : bipartitions2.keySet()) {
			if(!bipartitions1.containsKey(bv)) {
				bv.not();
				if(!bipartitions1.containsKey(bv)) {
					_fn++;
				}
				bv.not();
			}
		}

		// computed weighted forms
		if(t1.getRoot().getChildCount() == 2) {
			_num_iedges1 = t1.getNodeCount() - t1.getLeafCount() - 2;
		} else {
			_num_iedges1 = t1.getNodeCount() - t1.getLeafCount() - 1;
		}

		if(t2.getRoot().getChildCount() == 2) {
			_num_iedges2 = t2.getNodeCount() - t2.getLeafCount() - 2;
		} else {
			_num_iedges2 = t2.getNodeCount() - t2.getLeafCount() - 1;
		}

		if (_num_iedges1 == 0) {
			_wfp = 0;
		}
		else {
			_wfp = ((double) _fp) / ((double) _num_iedges1);
		}

		if (_num_iedges2 == 0) {
			_wfn = 0;
		}
		else {
			_wfn = ((double) _fn) / ((double) _num_iedges2);
		}

		// done
		return;
	}

	/**
	 * @return the number of internal edges in the first input tree.
	 */
	public int getNumInternalEdges1() {
		return _num_iedges1;
	}

	/**
	 * @return the number of internal edges in the second input tree.
	 */
	public int getNumInternalEdges2() {
		return _num_iedges2;
	}

	/**
	 * @return the number of false positives computed by the last call to {@link computeDifference}
	 */
	public int getFalsePositiveCount() {
		return _fp;
	}

	/**
	 * @return the number of false negatives computed by the last call to {@link computeDifference}
	 */
	public int getFalseNegativeCount() {
		return _fn;
	}

	public double getWeightedFalsePositive() {
		return _wfp;
	}

	public double getWeightedFalseNegative() {
		return _wfn;
	}

	public double getWeightedAverage() {
		return (_wfp + _wfn) / 2.0;
	}

	public double getUnweightedAverage() {
		return (_fp + _fn) / 2.0;
	}


}
