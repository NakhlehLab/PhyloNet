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

import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.recomp.WindowComparator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt.RiataHgt;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.*;
import java.util.Arrays;
import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 1/24/12
 * Time: 6:10 PM
 * To change this template use File | Settings | File Templates.
 */

public abstract class PAUPWindowComparator implements WindowComparator {

	// enumerations
	/**
	 * This enumeration is used to specify the standard tree distance
	 * measure that will be used by the window comparator
	 */
	public enum DistanceMeasure {
		SPR,
		RF
	}

	// constants
	protected static final String TAXA_PREFIX = "t";

	// inner classes
	public class PAUPTree {
		public MutableTree tree;
		public int score;

		public PAUPTree(MutableTree tree, int score) {
			this.tree = tree;
			this.score = score;
		}
	}

	// fields
	protected DistanceMeasure _measure;
	protected SymmetricDifference _sd;
	protected RiataHgt _riata;

	protected File _paup_path;
	protected int _num_iterations;
	protected int _pop_size;

	// constructors
	public PAUPWindowComparator(File paup_path, int pop_size, int num_iterations, DistanceMeasure dm) throws FileNotFoundException {
		_num_iterations = num_iterations;
		_pop_size = pop_size;

		// setup the distance measure
		_measure = dm;

		if(_measure == DistanceMeasure.RF) {
			_sd = new SymmetricDifference();
		} else {
			_riata = new RiataHgt();
		}

		// setup the paup path
		_paup_path = paup_path;
		if(_paup_path.exists() == false) {
			throw new FileNotFoundException("PAUP could not be found");
		}
	}

	// methods
	/**
	 * Compute the distance between two trees.
	 */
	protected double computeDistance(Tree t1, Tree t2) {
		if(_measure == DistanceMeasure.RF) {
			_sd.computeDifference(t1,t2,false);

			return _sd.getWeightedAverage();
		} else {
			 _riata.computeHgt(t1,t2);

			 // return the number of events
			 return _riata.getMinimumSolutionSize();
		}
	}

	protected PAUPTree[] computeMPTrees(String[] taxa_seqs)
		throws ParseException, IOException, InterruptedException {

		return computeMPTrees(_paup_path, taxa_seqs, _pop_size, _num_iterations);
	}

	/**
	 * Compute maximum parsimony trees for the specified set of sequences.
	 *
	 * @param paup_path is a valid path to the paup executable
	 * @param taxa_seqs is the set of sequences that will be used to compute the trees
	 * @param pop_size is the maximum number of trees to compute.
	 * @param num_iterations is the number of search iterations to use before stopping the heuristic search.
	 */
	protected PAUPTree[] computeMPTrees(File paup_path, String[] taxa_seqs, int pop_size, int num_iterations)
		throws ParseException, IOException, InterruptedException {

		// Setup PAUP Nexus file
		File nexus_file = File.createTempFile("nexus","paup");
		File tree_file = File.createTempFile("trees","out");
		File score_file = File.createTempFile("scores","out");

	//	CommandLineApp.printDebugMsg("\tNEXUS " + nexus_file.getAbsolutePath());
//		CommandLineApp.printDebugMsg("\tTREE " + tree_file.getAbsolutePath());
//		CommandLineApp.printDebugMsg("\tSCORE " + score_file.getAbsolutePath());

		writeMPNexusFile(nexus_file, tree_file, score_file, taxa_seqs, pop_size, num_iterations);

		// Run PAUP
		Process p = Runtime.getRuntime().exec(paup_path.getAbsolutePath() + " -n " + nexus_file.getAbsolutePath());
		p.waitFor();

		// Parse PAUP Resuls
		return parsePaupMPResults(tree_file, score_file, pop_size);
	}

	private PAUPTree[] parsePaupMPResults(File tree_file, File score_file, int pop_size) throws ParseException, IOException, FileNotFoundException {
		LineNumberReader treader = new LineNumberReader(new FileReader(tree_file));
		LineNumberReader sreader = new LineNumberReader(new FileReader(score_file));

		// skip the score header
		sreader.readLine();

		// read tree and score in pairs
		PAUPTree[] result = new PAUPTree[pop_size];
		int num_trees = 0;
		while(true) {

			String score_line = sreader.readLine();

			// if there are no more trees to read
			if(score_line == null) {
				break;
			}

			String[] score_info = score_line.trim().split("\\s+");

			result[num_trees] = new PAUPTree(new STITree<Object>(treader.readLine().trim()), Integer.parseInt(score_info[1]));
			num_trees++;
		}

		// "shorten" the array
		if(num_trees < pop_size) {
			PAUPTree[] tmp_result = new PAUPTree[num_trees];
			System.arraycopy(result,0,tmp_result,0,num_trees);
			result = tmp_result;
		}

		// sort the array
		Arrays.sort(result, new Comparator<PAUPTree>() {
            public int compare(PAUPTree h1, PAUPTree h2) {
                if (h1.score < h2.score) {
                    return -1;
                } else if (h1.score == h2.score) {
                    return 0;
                } else {
                    return 1;
                }
            }
        });

		return result;
	}

	private void writeMPNexusFile(File nexus_file, File tree_file, File score_file, String[] taxa_seqs, int pop_size, int num_iterations)
		throws FileNotFoundException {
		PrintStream ps = new PrintStream(new FileOutputStream(nexus_file));

		ps.println("#NEXUS");
		ps.println("begin paup;");
		ps.println("set Monitor=No MaxTrees=10000;");
		ps.println("end;");

		// write taxa names
		ps.println("begin taxa;");
		ps.println("dimensions ntax=" + taxa_seqs.length + ";");
		ps.println("taxlabels");
		for(int i = 0; i < taxa_seqs.length; i++) {
			ps.print(" " + TAXA_PREFIX + i);
		}
		ps.println(";");
		ps.println("end;");

		// write sequence data
		ps.println("begin characters;");
		ps.println("dimensions nchar=" + taxa_seqs[0].length() + ";");
		ps.println("format gap=- datatype=DNA;");
		ps.println("matrix");
		for(int i = 0; i < taxa_seqs.length; i++) {
			ps.println(TAXA_PREFIX + i + "\t" + taxa_seqs[i]);
		}
		ps.println(";");
		ps.println("end;");

		// write paup block
		ps.println("begin paup;");
		ps.println("set criterion=parsimony maxtrees=" + pop_size + " increase=no;");
		ps.println("hsearch start=stepwise addseq=random nreps=25 swap=tbr;");
		ps.println("filter best=yes;");
		ps.println("set maxtrees=100 increase=no;");
		ps.println("hsearch start=current swap=tbr hold=1 nreps=" + num_iterations + " nbest=" + pop_size + ";");
		ps.println("pscores / scorefile=" + score_file.getAbsolutePath() + ";");
		ps.println("savetrees file=" + tree_file.getAbsolutePath() + " format=phylip;");
		ps.println("quit;");
		ps.println("end;");

		ps.close();
	}

	protected int getEndOfTreeLevel(PAUPTree[] trees, int num_levels) {
		int curr_level = 0;
		int idx = 1;
		while(idx < trees.length) {
			if(trees[idx].score != trees[idx-1].score) {
				curr_level++;
			}

			if(curr_level == num_levels) {
				break;
			} else {
				idx++;
			}
		}

		return (idx - 1);
	}

	public boolean invertNormalizedScore() {
		return false;
	}

}
