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

package edu.rice.cs.bioinfo.programs.phylonet.ilp;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MaxClique;
import edu.rice.cs.bioinfo.programs.phylonet.algos.lca.SchieberVishkinLCA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/3/11
 * Time: 10:50 AM
 * To change this template use File | Settings | File Templates.
 */
public class IlpGenerator {

    private double _epsilon;	// To simulate strict inequalities in CPLEX.
	private double _scale;		// To scale in case the time is too small.
	private DecimalFormat _df;	// Decimal formater to print numbers in a particular format.
	private List<String> _taxa;	// Taxon names.
	private double _sfweight;	// Weight for sf (number of deep coalescences).
	private double _sgweight;	// Weight for sg (number of no coalescences).

	public IlpGenerator() {
		_epsilon = 1e-9;
		_scale = 1.0;
		_df = new DecimalFormat("######.#########");
		_sfweight = _sgweight = 1.0;
		_taxa = null;
	}

	public void setTaxa(List<String> t) {
		_taxa = t;
	}

	public static void genSpeciesTrees(String args[]) {
		if (args.length <= 0 || args[0].equals("-h")) {
			printUsage();
			return;
		}

		File gtfile = new File(args[0]);
		File stfile = new File(args[1]);

		List<Tree> geneTrees = new LinkedList<Tree>();

		try {
			BufferedReader br = new BufferedReader(new FileReader(gtfile));
			String line;

			while ((line = br.readLine()) != null) {
				if (line.length() > 0) {
					NewickReader nw = new NewickReader(new StringReader(line));
					geneTrees.add(nw.readTree());
				}
			}

			br.close();
		}
		catch (Exception e) {
			System.err.println(e.getMessage());
			return;
		}

		IlpGenerator builder = new IlpGenerator();

		try {
			builder.generateSpeciesTrees(geneTrees, stfile);
		} catch (IOException e) {
			// error while opening, writing to, or closing destination file
			System.err.println("I/O error while opening or writing to the destination file: " + e.getMessage());
			return;
		}

	}

	public static void genCplex(String args[]) {
		if (args.length <= 0 || args[0].equals("-h")) {
			printUsage();
			return;
		}

		String stfile = args[0];
		String gtfile = args[1];

		if (args.length > 2 && args.length < 4) {
			printUsage();
			return;
		}

		double sfw = Double.parseDouble(args[2]);
		double sgw = Double.parseDouble(args[3]);

		// Generate CPLEX inputs for each trees in the stfile.
		BufferedReader speciesReader, geneReader;

		try {
			speciesReader = new BufferedReader(new FileReader(stfile));
			geneReader = new BufferedReader(new FileReader(gtfile));
		}
		catch (Exception e) {
			System.err.println(e.getMessage());
			return;
		}

		List<STITree<Integer>> speciesTrees = new LinkedList<STITree<Integer>>();
		List<Tree> geneTrees = new LinkedList<Tree>();

		try {
			String line;
			while ((line = speciesReader.readLine()) != null) {
				NewickReader nw = new NewickReader(new StringReader(line));
				STITree<Integer> st = new STITree<Integer>(nw.readTree());

				speciesTrees.add(st);
			}

			while ((line = geneReader.readLine()) != null) {
				if (line.length() > 0) {
					NewickReader nw = new NewickReader(new StringReader(line));
					geneTrees.add(nw.readTree());
				}
			}

			speciesReader.close();
			geneReader.close();
		}
		catch (Exception e) {
			System.err.println(e.getMessage());
			return;
		}

		for (int i = 0; i < speciesTrees.size(); i++) {
			String cplexFileName = "input" + i;
			String variableFileName = "var" + i;
			String scriptFileName = "script" + i;

			IlpGenerator builder = new IlpGenerator();

			builder.setSfWeight(sfw);
			builder.setSgWeight(sgw);
			builder.generateCplexInput(speciesTrees.get(i), geneTrees, cplexFileName, variableFileName, scriptFileName);
		}
	}

	/**
	 * Set weight for sf.
	 *
	 * @param w: sf weight.
	 */
	public void setSfWeight(double w) {
		_sfweight = w;
	}

	/**
	 * Set weight for sg.
	 *
	 * @param w: sg weight.
	 */
	public void setSgWeight(double w) {
		_sgweight = w;
	}

	/**
	 * Generates a list of species trees by constructing from maximal sets of
	 * compatible clusters induced by gene trees, and writes this list to the
	 * supplied file.
	 *
	 * @param geneTrees		the list of gene trees from which we want to build
	 * 						species trees
	 * @param destination	the file to which the generated species trees
	 * 						should be written
	 *
	 * @throws IOException	if an I/O error is encountered while opening,
	 * 						writing to, or closing the destination file
	 */
	public void generateSpeciesTrees(List<Tree> geneTrees, File destination) throws IOException {

		// Generate all clusters in the gene trees.
		if (_taxa == null) {
			_taxa = new LinkedList<String>();
			for (TNode node : geneTrees.get(0).getNodes()) {
				if (node.isLeaf()) {
					_taxa.add(node.getName());
				}
			}
		}

		List<BitSet> clusters = computeAllClusters(geneTrees);

		// Create a graph of compatible clusters.
		double[][] compatibilityMatrix = new double[clusters.size()][clusters.size()];

		for (int i = 0; i < clusters.size() - 1; i++) {
			for (int j = i + 1; j< clusters.size(); j++) {
				BitSet bsi = clusters.get(i);
				BitSet bsj = clusters.get(j);

				if (areClustersCompatible(bsi, bsj)) {
					compatibilityMatrix[i][j] = 1;
				}
				else {
					compatibilityMatrix[i][j] = 0;
				}
			}
		}

		// Get all maximal cliques which correspond to maximal sets of compatible clusters.
		MaxClique mc = new MaxClique(compatibilityMatrix);
		List<int[]> nodeCliques = mc.calculateGroups(MaxClique.CLIQUES, 0);

		// Construct species trees from maximal cliques and write them to file.
		FileWriter fw = new FileWriter(destination);
		try {
			for (int[] nodeClique : nodeCliques) {
				List<BitSet> clique = new LinkedList<BitSet>();
				for (int nodeIndex : nodeClique) {
					clique.add(clusters.get(nodeIndex));
				}
				Tree t = buildTreeFromClusters(clique);
				fw.write(t.toNewick() + "\n"); // TODO Add tab before each tree to clarify start/stop when viewing with word wrap?
			}
		} finally {
			fw.close();
		}

	}

	/**
	 * This function generates CPLEX inputs according to the method described in the paper "Accurate and efficient
	 * tree reconstruction from genome-size multi-locus data under the coalescent".
	 *
	 * @param speciesTree: File containing the species tree
	 * @param geneTrees: File containing the gene trees. The function assumes that all gene trees in this file
	 * 					 satisfy the property that time(parent) > time(children).
	 * @param cplexFileName: File will contain the input to CPLEX
	 * @param variableFileName: File will contain variable mappings. This is valuable when one wants to parse CPLEX outputs.
	 * @param scriptFileName: File containing CPLEX commands to run CPLEX, and get the optimization solution.
	 */
	public void generateCplexInput(STITree<Integer> speciesTree, List<Tree> geneTrees, String cplexFileName, String variableFileName, String scriptFileName) {
		// 0. Get the set of leaves.
		if (_taxa == null) {
			_taxa = new LinkedList<String>();
			for (TNode node : speciesTree.getNodes()) {
				if (node.isLeaf()) {
					_taxa.add(node.getName());
				}
			}
		}

		// 1. Label internal nodes of the species tree, ie, give unique names to internal nodes.
		for (TNode node : speciesTree.getNodes()) {
			((STINode) node).setParentDistance(TNode.NO_DISTANCE);
		}
		Trees.autoLabelNodes((MutableTree) speciesTree);

		// 2. Label internal edges with numbers.
		FileWriter variableFile, scriptFile;
		FileWriter lpFile;

		try {
			variableFile = new FileWriter(variableFileName);
			scriptFile = new FileWriter(scriptFileName);
			lpFile = new FileWriter(cplexFileName + ".lp");
		}
		catch (IOException e) {
			System.err.println(e.getMessage());
			return;
		}

		// 3. Build the map of variables to names, and write them to the variable file.
		StringBuffer mapString = new StringBuffer();
		StringBuffer scriptString = new StringBuffer();
		StringBuffer sfString = new StringBuffer();			// The sum (fy + gy - 1).
		StringBuffer sgString = new StringBuffer();			// The sum gy.
		StringBuffer stString = new StringBuffer();			// The sum of branch lengths.

		mapString.append("TREE\n" + speciesTree.toString() + "\n");
		scriptString.append("read " + cplexFileName + ".lp\n");
		scriptString.append("mipopt\n");
		scriptString.append("display solution objective\n");
		scriptString.append("display solution variables sf\n");	// Value of the sum (fy + gy - 1).
		scriptString.append("display solution variables sg\n");	// Value of the sum gy.
		scriptString.append("display solution variables st\n");	// Value of the sum of branch lengths.

		Map<TNode, String> speciesVariableMap = new HashMap<TNode, String>();
		TreePostOrderTraversal traversal = new TreePostOrderTraversal(speciesTree.getRoot());
		int speciesCount = 0;	// Indicates the index of nodes in post order.

		// 3.1. Map for internal nodes of the species tree.
		mapString.append("SPECIES VARIABLES\n");
		for (TNode node : traversal) {
			if (!node.isLeaf()) {
				speciesCount++;
				speciesVariableMap.put(node, "t" + speciesCount);
				mapString.append("t" + speciesCount + ":" + node.getName() + "\n");
				scriptString.append("display solution variables t" + speciesCount + "\n");
			}
		}

		// 3.2. Create a special variable t for a new artifical root of the species tree.
		speciesCount++;
		String rootVar = new String("t" + speciesCount);
		mapString.append(rootVar + ":ROOT\n");
		scriptString.append("display solution variables " + rootVar + "\n");

		// 3.3. Build map for gene variables.
		List<Map<TNode, Integer>> geneVariableMaps = new LinkedList<Map<TNode, Integer>>();
		List<Integer> geneCounts = new LinkedList<Integer>();
		int gc = 0;

		for (Tree gt : geneTrees) {
			Map<TNode, Integer> map = new HashMap<TNode, Integer>();

			traversal = new TreePostOrderTraversal(gt.getRoot());
			for (TNode node : traversal) {
				boolean added = true;

				if (node.isLeaf() || node.getParentDistance() == 0) {
					added = false;
				}
				else {
					// Check if it has any child with equal time.
					for (TNode child : node.getChildren()) {
						if (child.getParentDistance() == node.getParentDistance()) {
							added = false;
							break;
						}
					}
				}

				if (added) {
					gc++;
					map.put(node, gc);
				}
			}

			geneCounts.add(gc);	// Used to count the number of internal nodes in each gene trees.
			geneVariableMaps.add(map);
		}

		// 4. Build the LP formula.
		double maxTime = 0.0;		// The maximum time value of the root.

		for (Tree gt : geneTrees) {
			if (maxTime < gt.getRoot().getParentDistance()) {
				maxTime = gt.getRoot().getParentDistance();
			}
		}

		maxTime *= _scale;
		_epsilon *= _scale;

		// 4.1. The objective function:
		traversal = new TreePostOrderTraversal(speciesTree.getRoot());
		StringBuffer objString = new StringBuffer("minimize " + _sfweight + " sf " + _sgweight + " sg " + " st\n");

		// 4.2. The string to get the value of st = sum of (tu - tv).
		for (TNode node : traversal) {
			if (!node.isLeaf()) {
				String tv = speciesVariableMap.get(node);
				double coef = node.getChildCount() - 1;

				stString.append(coef + " " + tv + " + ");
			}
		}
		stString.append(rootVar + " - " + _df.format(_scale) + " st = 0\n");

		// 4.3. The string to get the value of sf = sum of fy + gy - 1. and the string to get the value of sg.
		int totalGeneVar = 0;
		for (Map<TNode, Integer> map : geneVariableMaps) {
			for (Integer var : map.values()) {
				String tmp = "f" + var + " + g" + var;
				sfString.append(tmp + " + ");
				sgString.append("g" + var + " + ");
			}

			totalGeneVar += map.size();
		}
		sfString.delete(sfString.length() - 2, sfString.length());
		sfString.append("- sf = " + totalGeneVar + "\n");
		sgString.delete(sgString.length() - 2, sgString.length());
		sgString.append("- sg = 0\n");

		// 4.4. The temporal contraints on species tree's nodes:
		StringBuffer constraintString = new StringBuffer("subject to\n");

		for (TNode node : speciesVariableMap.keySet()) {
			if (!node.isLeaf()) {
				String tv = speciesVariableMap.get(node);
				String tu = node.isRoot() ? rootVar : speciesVariableMap.get(node.getParent());

				constraintString.append(tu + " - " + tv + " >= " + _df.format(_epsilon) + "\n");
			}
		}

		// 4.5. Temporal constraints for coalescent.
		mapString.append("GENE VARIABLES\n");

		int alphaCount = 1;
		StringBuffer boundString = new StringBuffer();

		for (Map<TNode, Integer> map : geneVariableMaps) {
			for (TNode node : map.keySet()) {
				// Get the lca for T(node).
				Set<TNode> geneSubleaves = getSubleaves(node);
				Set<TNode> speciesSubleaves = new HashSet<TNode>();

				for (TNode gl : geneSubleaves) {
					speciesSubleaves.add(speciesTree.getNode(gl.getName()));
				}

				SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(speciesTree);
				STINode<Integer> lca = (STINode<Integer>) lcaFinder.getLCA(speciesSubleaves);
				int c = (areCladesEqual(node, lca, _taxa) ? 1 : 0);

				mapString.append("f" + map.get(node) + ":" + c + "\n");
				scriptString.append("display solution variables f" + map.get(node) + "\n");

				// Add constraint for gy.
				double val = node.getParentDistance() * _scale - maxTime + _epsilon;
				constraintString.append(speciesVariableMap.get(lca) + " - " + _df.format(maxTime) + " g" + map.get(node) + " >= " + _df.format(val) + "\n");

				val = maxTime - node.getParentDistance() * _scale;
				constraintString.append(speciesVariableMap.get(lca) + " - " + _df.format(val) + " g" + map.get(node) + " <= " + _df.format(node.getParentDistance() * _scale) + "\n");

				// And, for each edge above the lca, add two constraints for coalescent.
				StringBuffer alphaSum = new StringBuffer();
				StringBuffer fConstraint = new StringBuffer("f" + map.get(node) + " - ");

				int m = 1;
				while (lca != null) {
					String tt = speciesVariableMap.get(lca);
					String a = "a" + (alphaCount++);
					val = node.getParentDistance() * _scale + maxTime;
					String temp = tt + " + " + _df.format(maxTime) + " " + a + " <= " + _df.format(val) + "\n";
					constraintString.append(temp);

					String ts;
					if (lca.isRoot()) {
						ts = rootVar;
						lca = lca.getParent();
					}
					else {
						lca = lca.getParent();
						ts = speciesVariableMap.get(lca);
					}
					val = node.getParentDistance() * _scale - maxTime;
					temp = ts + " - " + _df.format(maxTime) + " " + a + " >= " + _df.format(val) + "\n";
					constraintString.append(temp);

					alphaSum.append(a + " + ");
					fConstraint.append(m + " " + a + " - ");
					boundString.append(a + " ");

					m++;
				}

				// Add a constraints to enforce only one alpha will be one, while all other will be zero.
				alphaSum.append("g" + map.get(node) + " = 1\n");
				fConstraint.delete(fConstraint.length() - 3, fConstraint.length());
				fConstraint.append(" = 0\n");
				boundString.append("g" + map.get(node) + "\n");

				constraintString.append(alphaSum);
				constraintString.append(fConstraint);
			}
		}

		constraintString.append(stString);
		constraintString.append(sfString);
		constraintString.append(sgString);

		// 4.6. Constraint for integer values of alphas.
		constraintString.append("binaries\n");
		constraintString.append(boundString + "end\n");

		// 5. Save cplexString, mapString, and scriptString to files.
		try {
			lpFile.write(objString.toString());
			lpFile.write(constraintString.toString());
			lpFile.close();

			variableFile.write(mapString.toString());
			variableFile.close();

			scriptString.append("quit\n");
			scriptFile.write(scriptString.toString());
			scriptFile.close();
		}
		catch (Exception e) {
			System.err.println(e.getMessage());
			return;
		}
	}

	/**
	 * This function filters from the set of gene trees those whose node times satisfy: time(parent) > time(children).
	 *
	 * @param gtFile: File containing the original gene trees.
	 * @param validGtFile: File containing new "valid" gene trees (ie., those satisfy the description above).
	 */
	public void getValidGeneTrees(String gtFile, String validGtFile) {
		BufferedReader br;
		FileWriter fw;
		int invalidCount = 0;

		try {
			br = new BufferedReader(new FileReader(gtFile));
			fw = new FileWriter(validGtFile);
			String line;

			while ((line = br.readLine()) != null && line.length() > 0) {
				NewickReader nw = new NewickReader(new StringReader(line));
				Tree tr = nw.readTree();

				TreePostOrderTraversal traversal = new TreePostOrderTraversal(tr.getRoot());
				boolean valid = true;

				for (TNode node : traversal) {
					if (!node.isRoot()) {
						double tv = node.getParentDistance() * 1e9;
						double tu = node.getParent().getParentDistance() * 1e9;

						if (tu < tv) {
							invalidCount++;
							valid = false;
							break;
						}
					}
				}

				if (valid) {
					fw.write(line + "\n");
				}
			}

			System.out.println("Number of bad trees: " + invalidCount++);

			br.close();
			fw.close();
		}
		catch (Exception e) {
			System.err.println(e.getMessage());
			return;
		}
	}

	/**
	 * Convert the tree with branch lengths into the tree with time. The time of the leaves are assigned 0.
	 *
	 * @param tree: The tree with branch lengths.
	 */
	public void assignNodeTime(STITree<Double> tree) {
		TreePostOrderTraversal traversal = new TreePostOrderTraversal(tree.getRoot());

		for (TNode node : traversal) {
			if (node.isLeaf()) {
				((STINode<Double>) node).setData(0.0);
			}
			else {
				STINode<Double> child = ((STINode<Double>) node).getChildren().iterator().next();
				((STINode<Double>) node).setData(child.getData() + child.getParentDistance());
			}
		}
	}

	/**
	 * This function returns the set of leaves under a node.
	 *
	 * @param subroot: The node we want to compute the set of leaves under it.
	 * @return set of leaves under node subroot.
	 */
	private Set<TNode> getSubleaves(TNode subroot) {
		Set<TNode> subleaves = new HashSet<TNode>();
		TreePostOrderTraversal traversal = new TreePostOrderTraversal(subroot);

		for (TNode node : traversal) {
			if (node.isLeaf()) {
				subleaves.add(node);
			}
		}

		return subleaves;
	}

	/**
	 * This function checks if two clades are equal or not. (Two clades are equal if they have the same topology over
	 * the same set of taxa.) Two clades to be compared need not to be in the same tree.
	 *
	 * @param node1: Root of clade 1.
	 * @param node2: Root of clade 2.
	 * @return true if two clades are equal; false otherwise.
	 */
	private boolean areCladesEqual(TNode node1, TNode node2, List<String> taxa) {
		// 1. Get the set of all clusters under node1 and node2.
		Map<TNode, BitSet> clusterMap1 = new HashMap<TNode, BitSet>();
		Map<TNode, BitSet> clusterMap2 = new HashMap<TNode, BitSet>();

		TreePostOrderTraversal traversal = new TreePostOrderTraversal(node1);
		for (TNode child : traversal) {
			BitSet bs = new BitSet(taxa.size());

			if (child.isLeaf()) {
				int index = taxa.indexOf(child.getName());
				bs.set(index);

			}
			else {
				for (TNode grandChild : child.getChildren()) {
					bs.or(clusterMap1.get(grandChild));
				}
			}

			clusterMap1.put(child, bs);
		}

		traversal = new TreePostOrderTraversal(node2);
		for (TNode child : traversal) {
			BitSet bs = new BitSet(taxa.size());

			if (child.isLeaf()) {
				int index = taxa.indexOf(child.getName());
				bs.set(index);

			}
			else {
				for (TNode grandChild : child.getChildren()) {
					bs.or(clusterMap2.get(grandChild));
				}
			}

			clusterMap2.put(child, bs);
		}

		// 2. Check if both clades have the set of clusters.
		if (clusterMap1.get(node1).cardinality() != clusterMap2.get(node2).cardinality()) {
			return false;	// They have different set of labels, so they are different.
		}
		else if (clusterMap1.size() != clusterMap2.size()) {
			return false;	// They have different topology because they have different clusters.
		}
		else {
			for (BitSet bs1 : clusterMap1.values()) {
				if (!clusterMap2.containsValue(bs1)) {
					return false;
				}
			}
			for (BitSet bs2 : clusterMap2.values()) {
				if (!clusterMap1.containsValue(bs2)) {
					return false;
				}
			}

			return true;
		}
	}

	/**
	 * Print help message about how to use
	 */
	public static void printUsage() {
		System.out.println();
		System.out.println("To generate CPLEX inputs, type:");
		System.out.println("gencplex stfile gtfile");
		System.out.println("\tstfile: the file containing the species tree. The file can have more than one tree,");
		System.out.println("\t        in which case the program will generate CPLEX inputs for each species tree.");
		System.out.println("\tgtfile: the file containing the gene trees. All gene trees are assumed to satisfy the");
		System.out.println("\t        that time(parent) > time(children)");
		System.out.println();
		System.out.println("To generate species tree topologies, type:");
		System.out.println("genst gtfile stfile");
		System.out.println("\tgtfile: the file containing the gene trees.");
		System.out.println("\tstfile: the file to store generated species trees.");

/*		System.out.println();
		System.out.println("To convert the trees with branch lengths to trees with times, type:");
		System.out.println("convert file newfile");
		System.out.println("\tfile: the file that contains trees with branch lengths.");
		System.out.println("\tnewfile: the file you want to save converted trees (ie, those with times)");

		System.out.println();
		System.out.println("To filter out trees whose time values violate the property time(parent) > time(childrent), type:");
		System.out.println("filter file newfile");
		System.out.println("\tfile: the file that contains the original timed trees.");
		System.out.println("\tnewfile: the file you want to save valid trees (ie., those satisty the above property)");
*/
		System.out.println();
		System.out.println("To display this help message, type: -h");
	}

	/**
	 * Compute all clusters induced by a list of trees.
	 *
	 * @param geneTrees: List of gene trees.
	 * @return: List of clusters.
	 */
	public List<BitSet> computeAllClusters(List<Tree> geneTrees) {
		List<BitSet> allClusters = new LinkedList<BitSet>();

		for (Tree gt : geneTrees) {
			List<BitSet> cls = computeTreeClusters(gt);

			for (BitSet bs1 : cls) {
				boolean found = false;
				for (BitSet bs2 : allClusters) {
					if (areClustersEqual(bs1, bs2)) {
						found = true;
						break;
					}
				}

				if (!found) {
					allClusters.add(bs1);
				}
			}
		}

		return allClusters;
	}

	/**
	 * Computer clusters induced by a tree.
	 *
	 * @param tree: The tree we want to compute clusters.
	 * @return: A list of clusters induced by the tree.
	 */
	public List<BitSet> computeTreeClusters(Tree tree) {
		TreePostOrderTraversal traversal = new TreePostOrderTraversal(tree.getRoot());
		List<BitSet> clusters = new LinkedList<BitSet>();
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();

		for (TNode node : traversal) {
			BitSet bs = new BitSet(_taxa.size());

			if (node.isLeaf()) {
				int i = _taxa.indexOf(node.getName());
				bs.set(i);

				map.put(node, bs);
			}
			else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = map.get(child);
					bs.or(childCluster);
				}

				map.put(node, bs);
			}

			if (bs.cardinality() <= _taxa.size() - 1 && bs.cardinality() > 1) {
				clusters.add(bs);
			}
		}

		return clusters;
	}

	/**
	 * Check if two clusters are equal or not.
	 *
	 * @param bs1, bs2: Two clusters to be checked for equality.
	 */
	private boolean areClustersEqual(BitSet bs1, BitSet bs2) {
		for (int i = 0; i < _taxa.size(); i++) {
			if (bs1.get(i) != bs2.get(i)) {
				return false;
			}
		}

		return true;
	}

	/**
	 * Check the compatibility of two clusters. Two input clusters are assumed
	 * to be different. Return true if bs1 is compatible with bs2, false otherwise.
	 *
	 * @param bs1, bs2: Two clusters to be checked for compatibility.
	 */
	private boolean areClustersCompatible(BitSet bs1, BitSet bs2) {
		int oneCount = 0, zeroCount = 0;

		for (int i = 0; i < _taxa.size(); i++) {
			if (bs1.get(i)) {
				if (bs2.get(i)) {
					oneCount++;
				}
				else {
					zeroCount++;
				}
			}
		}

		if (oneCount == bs1.cardinality() || zeroCount == bs1.cardinality()) {
			return true;
		}

		oneCount = 0; zeroCount = 0;
		for (int i = 0; i < _taxa.size(); i++) {
			if (bs2.get(i)) {
				if (bs1.get(i)) {
					oneCount++;
				}
				else {
					zeroCount++;
				}
			}
		}

		if (oneCount == bs2.cardinality() || zeroCount == bs2.cardinality()) {
			return true;
		}

		return false;
	}

	private Tree buildTreeFromClusters(List<BitSet> clusterSet) {
		MutableTree tree = new STITree<Object>();

		// Create a big star tree.
		for (int i = 0; i < _taxa.size(); i++) {
			tree.getRoot().createChild(_taxa.get(i));
		}

		// Use clusters to refine the tree above.
		for (BitSet bs : clusterSet) {
			Set<TNode> leaves = new HashSet<TNode>();

			for (int i = 0; i < _taxa.size(); i++) {
				if (bs.get(i)) {
					TNode node = tree.getNode(_taxa.get(i));
					leaves.add(node);
				}
			}

			SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(tree);
			TNode lca = lcaFinder.getLCA(leaves);

			List<TNode> movedChildren = new LinkedList<TNode>();

			for (TNode child : lca.getChildren()) {
				BitSet childLeaves = new BitSet(_taxa.size());

				for (TNode leaf : getSubleaves(child)) {
					int index = _taxa.indexOf(leaf.getName());
					childLeaves.set(index);
				}

				//if (!areClustersCompatible(bs, childLeaves)) {
				int o = childLeaves.cardinality();
				childLeaves.and(bs);
				int n = childLeaves.cardinality();
				if (o == n) {
					movedChildren.add(child);
				}
			}

			STINode<Object> newChild = ((STINode<Object>) lca).createChild();

			while (!movedChildren.isEmpty()) {
				newChild.adoptChild((TMutableNode) movedChildren.get(0));
				movedChildren.remove(0);
			}
		}

		return tree;
	}
}

/**
 * An iterable class that allows us to visit the tree in the post order. To make the class more
 * general, we pass it a start node, so that it will visit all nodes under the start node in the
 * post order. If we want to visit all nodes in a tree, simply pass it the tree's root.
 */
class TreePostOrderTraversal implements Iterable<TNode> {
	private TNode _start_node;

	public TreePostOrderTraversal(TNode start) {
		_start_node = start;
	}

	public Iterator<TNode> iterator() {
		return new TreePostOrderIterator();
	}

	private class TreePostOrderIterator implements Iterator<TNode> {
		private Stack<TNode> _unvisited;		// Stores unvisited nodes.
		private Stack<Iterator<? extends TNode>> _iters;	// Stores pointers to children to be visited.

		public TreePostOrderIterator() {
			_unvisited = new Stack<TNode>();
			_iters = new Stack<Iterator<? extends TNode>>();

			if (_start_node != null) {
				_unvisited.push(_start_node);
				if (!_start_node.isLeaf()) {	// i.e., it has children.
					_iters.push(_start_node.getChildren().iterator());
				}
				else {
					_iters.push(null);
				}
			}
		}

		public boolean hasNext() {
			return !_unvisited.isEmpty();
		}

		public TNode next() {
			assert(!_unvisited.isEmpty());

			while (true) {
				Iterator<? extends TNode> it = _iters.peek();

				if (it != null && it.hasNext()) {	// Still has more children to visit.
					TNode child = it.next();

					_unvisited.push(child);
					if (!child.isLeaf()) {
						_iters.push(child.getChildren().iterator());
					}
					else {
						_iters.push(null);
					}
				}
				else {
					break;	// Found next node in post order traversal.
				}
			}

			_iters.pop();
			return _unvisited.pop();
		}

		public void remove() {
			System.err.println("This method is currently not supported.");
			return;
		}
	}
}
