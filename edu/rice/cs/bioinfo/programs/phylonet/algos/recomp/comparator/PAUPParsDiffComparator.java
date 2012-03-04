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

import edu.rice.cs.bioinfo.programs.phylonet.algos.fitchpars.ParsimonyCalculator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 1/25/12
 * Time: 6:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class PAUPParsDiffComparator extends PAUPWindowComparator {

	// fields
	private int _num_levels;

	private Hashtable<Integer,PAUPTree[]> _past_trees = new Hashtable<Integer,PAUPTree[]>();

	/**
	 * Construct a comparator that uses the top <code>num_levels</code> of the parsimony
	 * trees for a window.
	 *
	 * @param paup_path is the path to the PAUP executable on this system.
	 * @param num_levels the number of top levels to incorporate in the score computation.
	 * @param num_iterations the number of iterations the MP heuristic should use to find
	 * the optimal trees.
	 * @param pop_size the number of trees used in the MP search to find the optimal trees.
	 *
	 * @throws FileNotFoundException if paup_path does not specify a valid paup file
	 */
	public PAUPParsDiffComparator(File paup_path, int pop_size, int num_iterations, int num_levels)
		throws FileNotFoundException {

		super(paup_path, pop_size, num_iterations,DistanceMeasure.RF);

		_num_levels = num_levels;
	}

	public double compare(int window1_coord, int window2_coord, char[][] window1, char[][] window2) {

		String[] window1_strs = new String[window1.length];
		String[] window2_strs = new String[window2.length];

		for(int i = 0; i < window1.length; i++) {
			window1_strs[i] = new String(window1[i]);
			window2_strs[i] = new String(window2[i]);
		}

		// compute the parsimony trees for window1
		PAUPTree[] left_trees = null;

		// if we don't have the left window tree, make it - this should just happen at first.
		if(!_past_trees.containsKey(window1_coord)) {
			try {
		//		CommandLineApp.printDebugMsg("\tBuilding left trees");
				left_trees = computeMPTrees(window1_strs);
			//	CommandLineApp.printDebugMsg("\tBuilt left trees");
			} catch(IOException ioe) {
				System.err.println("ERROR: PAUPParsDiffComparator encountered I/O error");
				return Double.NaN;
			} catch(InterruptedException ie) {
				System.err.println("ERROR: PAUPParsDiffComparator interrupted");
				return Double.NaN;
			} catch(ParseException pe) {
				System.err.println("ERROR: PAUPParsDiffComparator encountered illegal PAUP tree");
				return Double.NaN;
			//} catch(ParseException pe) {
				//System.err.println("ERROR: PAUP returned an invalid tree");
			//	return Double.NaN;
			}
		} else {
			left_trees = _past_trees.remove(window1_coord);
			//CommandLineApp.printDebugMsg("\tRetrieved left trees");
		}



		// compute the parsimony trees for window2
		PAUPTree[] right_trees = null;
		try {
		//	CommandLineApp.printDebugMsg("\tBuilding right trees");
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
	//	} catch(ParseException pe) {
	//		System.err.println("ERROR: PAUP returned an invalid tree");
		//	return Double.NaN;
		}

	//	CommandLineApp.printDebugMsg("\tBuilt right trees");

		// generate names list
		String[] names = new String[window1.length];
		for(int i = 0; i < names.length; i++) {
			names[i] = TAXA_PREFIX + i;
		}

		// get the min pscore for the 2nd window
		int min_window2_pscore = right_trees[0].score;

		//int right_last_tidx = getEndOfTreeLevel(trees, _num_levels);
		int left_last_tidx = getEndOfTreeLevel(left_trees, _num_levels);

		// get the min pscore for the 1st window
	//	CommandLineApp.printDebugMsg("\tComputing parsimony scores");

		ParsimonyCalculator pcalc = new ParsimonyCalculator();
		int min_window1_pscore = Integer.MAX_VALUE;

		for(int idx = 0; idx <= left_last_tidx; idx++) {
			int pscore = pcalc.computeParsimony(left_trees[idx].tree, names, window2_strs);

			if(pscore < min_window1_pscore) {
				min_window1_pscore = pscore;
			}
		}

		// store the right trees for later use
		_past_trees.put(window2_coord, right_trees);

		// done
	//	CommandLineApp.printDebugMsg("\tDone");
		return Math.abs(min_window2_pscore - min_window1_pscore);
	}

	@Override
	public String toString() {
		return "Maximum Parsimony";
	}
}
