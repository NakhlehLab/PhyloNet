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

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Hashtable;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 1/25/12
 * Time: 6:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class PAUPAvgMaxComparator extends PAUPWindowComparator {

    // fields
	protected int _num_levels;

	protected Hashtable<Integer,PAUPTree[]> _past_trees = new Hashtable<Integer,PAUPTree[]>();

	// constructors
	public PAUPAvgMaxComparator(File paup_path, int pop_size, int num_iterations, int num_levels, PAUPWindowComparator.DistanceMeasure dm) throws FileNotFoundException {
		super(paup_path, pop_size, num_iterations, dm);

		_num_levels = num_levels;
	}
	// methods
	public double compare(int window1_coord, int window2_coord, char[][] window1, char[][] window2) {

		PAUPTree[] left_trees = _past_trees.remove(window1_coord);
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
	//	} catch(phylonet.tree.io.ParseException pe) {
//			System.err.println("ERROR: PAUP returned an invalid tree");
//			return Double.NaN;
		}

		// get the top levels
		int last_left_tree = getEndOfTreeLevel(left_trees, _num_levels);
		int last_right_tree = getEndOfTreeLevel(right_trees, _num_levels);

		// compute dist_max for all left trees
		double[] left_max_dists = new double[last_left_tree+1];
		double[] right_max_dists = new double[last_right_tree+1];

		Arrays.fill(left_max_dists, 0.0);
		Arrays.fill(right_max_dists, 0.0);

		for(int i = 0; i <= last_left_tree; i++) {
			for(int j = 0; j <= last_right_tree; j++) {
				double dist = computeDistance(left_trees[i].tree, right_trees[j].tree);

				if(dist > left_max_dists[i]) {
					left_max_dists[i] = dist;
				}

				if(dist > right_max_dists[j]) {
					right_max_dists[j] = dist;
				}
			}
		}

		double left_avg = 0;
		for(int i = 0; i < left_max_dists.length; i++) {
			left_avg += left_max_dists[i];
		}
		left_avg = left_avg / ((double) left_max_dists.length);

		double right_avg = 0;
		for(int i = 0; i < right_max_dists.length; i++) {
			right_avg += right_max_dists[i];
		}
		right_avg = right_avg / ((double) right_max_dists.length);

		// store the right tree
		_past_trees.put(window2_coord, right_trees);

		return (left_avg + right_avg);
	}

	@Override
	public String toString() {
		return "Average Maximum Distance";
	}
}
