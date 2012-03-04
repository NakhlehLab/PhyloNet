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

package edu.rice.cs.bioinfo.programs.phylonet.algos.recomp.comparator;

import edu.rice.cs.bioinfo.programs.phylonet.algos.consensus.TreeConsensusCalculator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 1/25/12
 * Time: 6:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class PAUPIntersectionComparator extends PAUPWindowComparator {

    // fields
	protected int _num_levels;

	protected Hashtable<Integer,PAUPWindowComparator.PAUPTree[]> _past_trees = new Hashtable<Integer,PAUPWindowComparator.PAUPTree[]>();

	// constructors
	public PAUPIntersectionComparator(File paup_path, int pop_size, int num_iterations, int num_levels, PAUPWindowComparator.DistanceMeasure dm) throws FileNotFoundException {
		super(paup_path, pop_size, num_iterations, dm);

		_num_levels = num_levels;
	}

	public double compare(int window1_coord, int window2_coord, char[][] window1, char[][] window2) {

		PAUPWindowComparator.PAUPTree[] left_trees = _past_trees.remove(window1_coord);
		PAUPWindowComparator.PAUPTree[] right_trees;

		String[] window1_strs = new String[window1.length];
		String[] window2_strs = new String[window2.length];

		for(int i = 0; i < window1.length; i++) {
			window1_strs[i] = new String(window1[i]);
			window2_strs[i] = new String(window2[i]);
		}

		// if we didn't get the left tree, compute it
		if(left_trees == null) {
			try {
				left_trees = computeMPTrees(window1_strs);
			} catch(IOException ioe) {
				System.err.println("ERROR: PAUPParsDiffComparator encountered I/O error");
				return Double.NaN;
			} catch(InterruptedException ie) {
				System.err.println("ERROR: PAUPParsDiffComparator interrupted");
				return Double.NaN;
			} catch(ParseException pe) {
				System.err.println("ERROR: PAUPParsDiffComparator encountered illegal PAUP tree");
				return Double.NaN;
		//	} catch(phylonet.tree.io.ParseException pe) {
		//		System.err.println("ERROR: PAUP returned an invalid tree");
		//		return Double.NaN;
			}
		}

		// compute the parsimony trees for window2
		try {
			right_trees = computeMPTrees(window2_strs);
		} catch(IOException ioe) {
			System.err.println("ERROR: PAUPParsDiffComparator encountered I/O error");
			return Double.NaN;
		} catch(InterruptedException ie) {
			System.err.println("ERROR: PAUPParsDiffComparator interrupted");
			return Double.NaN;
		} catch(ParseException pe) {
			System.err.println("ERROR: PAUPParsDiffComparator encountered illegal PAUP tree");
			return Double.NaN;
		//} catch(phylonet.tree.io.ParseException pe) {
	//		System.err.println("ERROR: PAUP returned an invalid tree");
	//		return Double.NaN;
		}

		// get the top levels
		int last_left_tree = getEndOfTreeLevel(left_trees, _num_levels);
		int last_right_tree = getEndOfTreeLevel(right_trees, _num_levels);

		// compute consensus trees
		TreeConsensusCalculator tcc = new TreeConsensusCalculator();
		Set<Tree> trees = new HashSet<Tree>();

		for(int i = 0; i <= last_left_tree; i++) {
			trees.add(left_trees[i].tree);
		}

		Tree left_consensus = tcc.computeUnrootedConsensus(trees, 1.0);

		trees.clear();
		for(int i = 0; i <= last_right_tree; i++) {
			trees.add(right_trees[i].tree);
		}

		Tree right_consensus = tcc.computeUnrootedConsensus(trees, 1.0);

		// compute radii
		double left_radius = computeRadius(left_consensus, left_trees, last_left_tree);
		double right_radius = computeRadius(right_consensus, right_trees, last_right_tree);

		// count trees closer
		int num_close_rtrees = countCloseTrees(left_consensus, left_radius, right_trees, last_right_tree);
		int num_close_ltrees = countCloseTrees(right_consensus, right_radius, left_trees, last_left_tree);

		double dist1 = ((double) num_close_rtrees) / ((double) (last_right_tree+1));
		double dist2 = ((double) num_close_ltrees) / ((double) (last_left_tree+1));

		// save the right trees for re-use later
		_past_trees.put(window2_coord, right_trees);

		// done
		return dist1 + dist2;
	}

	/**
	 * Count the number of trees that are within <code>radius</code> of the consensus tree.
	 */
	protected int countCloseTrees(Tree consensus, double radius, PAUPTree[] trees, int last_tree) {

		int num_trees = 0;

		for(int i = 0; i <= last_tree; i++) {
			if(computeDistance(consensus, trees[i].tree) <= radius) {
				num_trees++;
			}
		}

		return num_trees;
	}

	/**
	 * Compute the largest distance that any tree is from the consensus tree.  This distance is called
	 * the radius.
	 */
	protected double computeRadius(Tree consensus, PAUPTree[] trees, int last_tree) {

		double radius = 0;

		for(int i = 0; i <= last_tree; i++) {
			double dist = computeDistance(consensus, trees[i].tree);

			if(dist > radius) {
				radius = dist;
			}
		}

		return radius;
	}

	public boolean invertNormalizedScore() {
		return true;
	}

	@Override
	public String toString() {
		return "Percent Intersection";
	}
}
