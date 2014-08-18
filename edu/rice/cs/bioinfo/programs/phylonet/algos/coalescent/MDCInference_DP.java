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

package edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MaxClique;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Collapse;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.PostTraversal;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/9/11
 * Time: 3:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class MDCInference_DP
{
	/**
	 * Infers the species tree from the given list of rooted gene trees with single allele.
	 *
	 * @param 	trees	the list of gene trees from which the species tree is to be inferred
	 * @param	explore	false if the result is the only one optimal tree
	 *                  true if the result is a set of trees
	 * @param   proportion	  specified when explore is set to true
	 *          			  specifies the range of the number of extra lineages that the set of trees have
	 * @param   exhaust	  false if clusters are induced from gene trees
	 * 					  true if clusters all all possible clusters
	 * @return	tree(s) inferred from the supplied gene trees with the total number of extra lineages
	 */
	public List<Solution> inferSpeciesTree(List<MutableTuple<Tree,Double>> trees, boolean explore, double proportion, boolean exhaust, boolean unresolved, double time) {

		if (trees == null || trees.size() == 0) {
			throw new IllegalArgumentException("empty or null list of trees");
		}


		//Collapse.CollapseDescriptor cd = doCollapse(trees);
		//System.out.println("after collapse");

		List<String> taxalist = new ArrayList<String>();
        for(MutableTuple<Tree,Double> tr: trees){
			for (TNode node : tr.Item1.postTraverse()) {
				if (node.isLeaf() && !taxalist.contains(node.getName())) {
					taxalist.add(node.getName());
				}
			}
		}

		String[] taxa = new String[taxalist.size()];

		int index = 0;
		for(String taxon: taxalist){
			taxa[index++] = taxon;
		}


		//String[] taxa = trees.get(0).getLeaves();
		// Find the tree with minimum score over all taxa.
		Map<Integer, List<Vertex>> clusters = new HashMap<Integer, List<Vertex>>();

		double maxEL;
		if(!exhaust){
			maxEL = computeTreeClusters(trees, taxa, clusters);
		}
		else{
			maxEL = computeAllClusters(trees, taxa, clusters);
		}
		List<Solution> solutions;
		//System.out.println("after collapse cluster");
		if(explore){
			solutions = findTreesByClique(clusters,taxa,proportion);
		}
		else{
			solutions = findTreesByDP(clusters,taxa,maxEL);
		}
		//System.out.println("after dp");
		if(!unresolved){
			time = time * 60;
			for(Solution sol: solutions){
				if(!Trees.isBinary(sol._st)){
					sol._totalCoals = tryBinaryResolutions(sol._st, time, taxa, trees, null) + sol._totalCoals;
				}
			}
		}

		//restoreCollapse(solutions, trees, cd);
		return solutions;
	}


	/**
	 * Infers the species tree from the given list of rooted gene trees with multiple alleles.
	 *
	 * @param 	trees	the list of gene trees from which the species tree is to be inferred
	 * @param	taxonMap	Maps gene tree taxa to species tree taxa.
	 * @param	explore	false if the result is the only one optimal tree
	 *                  true if the result is a set of trees
	 * @param   proportion	  specified when explore is set to true
	 *          			  specifies the range of the number of extra lineages that the set of trees have
	 * @param   exhaust	  false if clusters are induced from gene trees
	 * 					  true if clusters all all possible clusters
	 * @return	tree(s) inferred from the supplied gene trees with the total number of extra lineages
	 */

	//TODO cd is not correct for taxonmap
	public List<Solution> inferSpeciesTree(List<MutableTuple<Tree,Double>> trees, Map<String, String> taxonMap, boolean explore, double proportion, boolean exhaust, boolean unresolved, double time) {
		if (trees == null || trees.size() == 0) {
			System.err.println("Empty list of trees. The function returns a null tree.");
			return null;
		}

		String error = Trees.checkMapping(trees, taxonMap);
		if(error!=null){
			throw new RuntimeException("Gene trees have leaf named " + error + " that hasn't been defined in the mapping file");
		}

		List<String> temp1 = new LinkedList<String>();
		List<String> temp2 = new LinkedList<String>();
		for (String s : taxonMap.keySet()) {
			temp1.add(s);	// Gene tree taxa.
			if (!temp2.contains(taxonMap.get(s))) {
				temp2.add(taxonMap.get(s));	// Species tree taxa.
			}
		}

		String gtTaxa[] = new String[temp1.size()];
		String stTaxa[] = new String[temp2.size()];

		for (int i = 0; i < gtTaxa.length; i++) {
			gtTaxa[i] = temp1.get(i);
		}
		for (int i = 0; i < stTaxa.length; i++) {
			stTaxa[i] = temp2.get(i);
		}

		// Find the tree with the minimum score.
		Map<Integer, List<Vertex>> clusters = new HashMap<Integer, List<Vertex>>();
		double maxEL;
		if(!exhaust){
			maxEL = computeTreeClusters(trees, stTaxa, gtTaxa, taxonMap, clusters);
		}
		else{
			maxEL = computeAllClusters(trees, stTaxa, taxonMap, clusters);
		}
		List<Solution> solutions;
		if(explore){
			solutions = findTreesByClique(clusters,stTaxa,proportion);
		}
		else{
			solutions = findTreesByDP(clusters,stTaxa, maxEL);
		}

		if(!unresolved){
			time = time * 60;
			for(Solution sol: solutions){
				if(!Trees.isBinary(sol._st)){
					sol._totalCoals = tryBinaryResolutions(sol._st, time, stTaxa, trees, taxonMap) + sol._totalCoals;
				}
			}
		}

		return solutions;
	}


	private Collapse.CollapseDescriptor doCollapse(List<Tree> trees){
		Collapse.CollapseDescriptor cd = Collapse.collapse(trees);
		return cd;
	}

	private void restoreCollapse(List<Solution> sols, List<Tree> gts, Collapse.CollapseDescriptor cd){
		for(Solution sol: sols){
			Tree tr = sol._st;
			Collapse.expand(cd, (MutableTree)tr);
			for (TNode node : tr.postTraverse()) {
				if(((STINode<Integer>)node).getData()==null){
					((STINode<Integer>)node).setData(0);
				}
			}
		}
        for(Tree gt: gts){
            Collapse.expand(cd, (MutableTree)gt);
        }
	}

	/**
	 * Find the species tree by dynamic programming
	 *
	 * @param 	clusters	clusters with extra lineages
	 * @param	stTaxa	  taxa expression
	 * @param   maxEL	  the maximal number of extra lineage
	 * @return	the tree inferred from the supplied gene trees with the total number of extra lineages
	 *
	 */
	private List<Solution> findTreesByDP(Map<Integer, List<Vertex>> clusters, String[] stTaxa, double maxEL){
		List<Solution> solutions = new ArrayList<Solution>();

		Vertex all = clusters.get(stTaxa.length).get(0);
		computeMinCost(clusters,all,maxEL);

		// Build the minimum tree.
		List<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		List<Double> coals = new LinkedList<Double>();
		Stack<Vertex> minVertices = new Stack<Vertex>();
		if (all._min_rc != null) {
			minVertices.push(all._min_rc);
		}
		if (all._min_lc != null) {
			minVertices.push(all._min_lc);
		}
		if (all._subcl != null) {
			for(Vertex v: all._subcl){
				minVertices.push(v);
			}
		}

		while (!minVertices.isEmpty()) {
			Vertex pe = minVertices.pop();

			minClusters.add(pe._cluster);
			coals.add(pe._el_num);
			if (pe._min_rc != null) {
				minVertices.push(pe._min_rc);
			}
			if (pe._min_lc != null) {
				minVertices.push(pe._min_lc);
			}
			if (pe._subcl != null){
				for(Vertex v: pe._subcl){
					minVertices.push(v);
				}
			}
		}
		Solution sol = new Solution();
		if (minClusters == null || minClusters.isEmpty()) {
			MutableTree tr = new STITree<Object>();
			for (String s : stTaxa) {
				tr.getRoot().createChild(s);
			}
			sol._st = tr;
		}
		else{
			sol._st = Trees.buildTreeFromClusters(minClusters);
		}

		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		for (TNode node : sol._st.postTraverse()) {
			BitSet bs = new BitSet();
			if (node.isLeaf()) {
				for (int i = 0; i < stTaxa.length; i++) {
					if (node.getName().equals(stTaxa[i])) {
						bs.set(i);
						break;
					}
				}
				map.put(node, bs);
			}
			else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = map.get(child);
					bs.or(childCluster);
				}
				map.put(node, bs);
			}
			STITreeCluster c = new STITreeCluster(stTaxa);
			c.setCluster(bs);
			if(c.getClusterSize()==stTaxa.length){
				((STINode<Double>)node).setData(0.0);
			}
			else{
				int pos = minClusters.indexOf(c);
				((STINode<Double>)node).setData(coals.get(pos));
			}
		}

		sol._totalCoals = all._min_cost;
		solutions.add(sol);

		return solutions;
	}


	/**
	 * Find the species tree by finding the maximal cliques
	 *
	 * @param 	cmap	clusters with extra lineages
	 * @param	stTaxa	  taxa expression
	 * @param   proportion	the range of the number of extra lineages that the set of trees have
	 * @return	the tree inferred from the supplied gene trees with the total number of extra lineages
	 *
	 */
	private List<Solution> findTreesByClique(Map<Integer, List<Vertex>> cmap, String[] stTaxa, double proportion){
		List<Solution> solutions = new LinkedList<Solution>();

		List<STITreeCluster<Double>> clusters = new ArrayList<STITreeCluster<Double>>();
		int addEL = 0;

		for(Map.Entry<Integer, List<Vertex>> entry: cmap.entrySet()){
			if(entry.getKey()==1){
				List<Vertex> l = entry.getValue();
				for(Vertex v: l){
					addEL += v._el_num;
				}
			}
			else if(entry.getKey() < stTaxa.length){
				List<Vertex> l = entry.getValue();
				for(Vertex v: l){
					STITreeCluster<Double> c = new STITreeCluster<Double>(v._cluster);
					c.setData(v._el_num);
					clusters.add(c);
				}
			}
		}

		double[][] compatibilityMatrix = new double[clusters.size()][clusters.size()];
		for (int i = 0; i < clusters.size() - 1; i++) {
			STITreeCluster cl1 = clusters.get(i);
			for (int j = i + 1; j< clusters.size(); j++) {
				STITreeCluster cl2 = clusters.get(j);

				if (cl1.isCompatible(cl2)) {
					compatibilityMatrix[i][j] = 1;
				}
				else {
					compatibilityMatrix[i][j] = 0;
				}
			}
		}

		// Get all maximal cliques and score them
		MaxClique mc = new MaxClique(compatibilityMatrix);
		List<int[]> nodeCliques = mc.calculateGroups(MaxClique.CLIQUES, 0);
		List<Solution> maxCliques = new LinkedList<Solution>();

		for(int max = clusters.get(0).getTaxa().length-2 ; max>1; max--){
			for(int[] nodes: nodeCliques){
				if(nodes.length != max){
					continue;
				}
				int sum = 0;
				for(int id: nodes){
					sum += clusters.get(id).getData();
				}
				Solution s = new Solution();
				s._clusterIDs = nodes;
				s._totalCoals = sum + addEL;
				maxCliques.add(s);
			}
			if(maxCliques.size()>0){
				break;
			}
		}

		//order the maximal cliques
		for(int i=1;i<maxCliques.size();i++){
			Solution s1 = maxCliques.get(i);
			for(int j=0;j<i;j++){
				Solution s2 = maxCliques.get(j);
				if(s1._totalCoals < s2._totalCoals){
					maxCliques.remove(s1);
					maxCliques.add(j,s1);
					break;
				}
			}
		}

		//get the result
		double minCoal = maxCliques.get(0)._totalCoals;
		double maxCoal = (int)(( 1 + proportion / 100 ) * minCoal);
		for(Solution s: maxCliques){
			if(s._totalCoals <= maxCoal){
				List<STITreeCluster> minClusters = new ArrayList<STITreeCluster>();
				List<Double> coals = new ArrayList<Double>();
				for(int id: s._clusterIDs){
					STITreeCluster c = clusters.get(id);
					minClusters.add(c);
					coals.add(((STITreeCluster<Double>)c).getData());
				}
				for(Vertex v: cmap.get(1)){
					minClusters.add(v._cluster);
					coals.add(v._el_num);
				}
				Tree st = Trees.buildTreeFromClusters(minClusters);
				Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
				for (TNode node : st.postTraverse()) {
					BitSet bs = new BitSet();
					if (node.isLeaf()) {
						for (int i = 0; i < stTaxa.length; i++) {
							if (node.getName().equals(stTaxa[i])) {
								bs.set(i);
								break;
							}
						}
						map.put(node, bs);
					}
					else {
						for (TNode child : node.getChildren()) {
							BitSet childCluster = map.get(child);
							bs.or(childCluster);
						}
						map.put(node, bs);
					}
					STITreeCluster c = new STITreeCluster(stTaxa);
					c.setCluster(bs);
					if(c.getClusterSize() == stTaxa.length){
						((STINode<Double>)node).setData(0.0);
					}
					else{
						int pos = minClusters.indexOf(c);
						((STINode<Double>)node).setData(coals.get(pos));
					}
				}
				s._st = st;
				solutions.add(s);
			}
			else{
				break;
			}
		}
		return solutions;
	}


	/**
	 * This function infers the species tree using all possible clusters.
	 *
	 * @param 	trees	the list of gene trees
	 * @param	taxonMap	Maps gene tree taxa to species tree taxa.
	 * @return
	 */
	/*
	public Tree exhaustInferSpeciesTree(List<Tree> trees, Map<String, String> taxonMap) {
		Map<Integer, List<Vertex>> clusters = new HashMap<Integer, List<Vertex>>();
		if (trees == null || trees.size() == 0) {
			System.err.println("Empty list of trees. The function returns a null tree.");
			return null;
		}

		List<String> temp1 = new LinkedList<String>();
		List<String> temp2 = new LinkedList<String>();
		for (String s : taxonMap.keySet()) {
			temp1.add(s);	// Gene tree taxa.
			if (!temp2.contains(taxonMap.get(s))) {
				temp2.add(taxonMap.get(s));	// Species tree taxa.
			}
		}

		String gtTaxa[] = new String[temp1.size()];
		String stTaxa[] = new String[temp2.size()];

		for (int i = 0; i < gtTaxa.length; i++) {
			gtTaxa[i] = temp1.get(i);
		}
		for (int i = 0; i < stTaxa.length; i++) {
			stTaxa[i] = temp2.get(i);
		}

		// Find the tree with the minimum score.
		// computeTreeClusters(trees, stTaxa, gtTaxa, taxonMap);
		Integer maxEL = 0;
		generateClusters(trees, stTaxa, taxonMap, maxEL);

		Vertex all = clusters.get(stTaxa.length).get(0);
		computeMinCost(clusters,all, maxEL);

		// Build the minimum tree.
		List<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		Stack<Vertex> minVertices = new Stack<Vertex>();
		if (all._min_rc != null) {
			minVertices.push(all._min_rc);
		}
		if (all._min_lc != null) {
			minVertices.push(all._min_lc);
		}

		while (!minVertices.isEmpty()) {
			Vertex pe = minVertices.pop();

			minClusters.add(pe._cluster);
			if (pe._min_rc != null) {
				minVertices.push(pe._min_rc);
			}
			if (pe._min_lc != null) {
				minVertices.push(pe._min_lc);
			}
		}

		return Trees.buildTreeFromClusters(minClusters);
	}

	*/

	/**
	 * This function is written to compute clusters induced by gene trees when there are multiple
	 * gene copies.
	 *
	 * @param 	trees	the list of gene trees
	 * @param	stTaxa	taxa in species tree
	 * @param	gtTaxa	taxa in gene trees
	 * @param	taxonMap	Maps gene tree taxa to species tree taxa.
	 * @param	clusters	clusters with the number of extra lineages
	 * @return  the maximal number of extra lineages
	 */
	private double computeTreeClusters(List<MutableTuple<Tree,Double>> trees, String stTaxa[], String gtTaxa[], Map<String, String> taxonMap, Map<Integer, List<Vertex>> clusters) {
		double maxEL = 0;
        int total = 0;
        for(MutableTuple<Tree,Double> tr: trees){
            boolean contain = false;
			for (STITreeCluster tc : tr.Item1.getClusters(gtTaxa, false)) {
				// Create clusters that do not contain multiple copies of a taxon.
				STITreeCluster stCluster = new STITreeCluster(stTaxa);
				for (String s : tc.getClusterLeaves()) {
					stCluster.addLeaf(taxonMap.get(s));
				}

                if(stCluster.containsLeaf("S") && stCluster.containsLeaf("M") && stCluster.getClusterSize()==2){
                    System.out.println(tc);
                    contain = true;
                }

				// Added it to the list.
				int csize = stCluster.getClusterSize();

				if (clusters.containsKey(csize)) {
					List<Vertex> l = clusters.get(csize);
					boolean found = false;

					for (Vertex v : l) {
						if (v._cluster.equals(stCluster)) {
							found = true;
							break;
						}
					}

					if (!found) {
						Vertex nv = new Vertex();
						nv._cluster = stCluster;
						nv._el_num = DeepCoalescencesCounter.getClusterCoalNum(trees, stCluster, taxonMap, true);
						if(nv._el_num > maxEL){
							maxEL = nv._el_num;
						}
						nv._min_cost = -1;

						l.add(nv);
					}
				}
				else if (csize > 1) {	// Clusters of size 1 will be added later.
					List<Vertex> l = new LinkedList<Vertex>();
					Vertex v = new Vertex();

					v._cluster = stCluster;
					v._el_num = DeepCoalescencesCounter.getClusterCoalNum(trees, stCluster, taxonMap, true);
					if(v._el_num > maxEL){
						maxEL = v._el_num;
					}
					v._min_cost = -1;

					l.add(v);
					clusters.put(csize, l);
				}
			}
            if(contain){
                total += tr.Item2;
            }
		}

		// Add the cluster containing all taxa to the map.
		STITreeCluster all = new STITreeCluster(stTaxa);
		for (String t : stTaxa) {
			all.addLeaf(t);
		}

		Vertex v = new Vertex();
		v._cluster = all;
		v._el_num = 0;
		v._min_cost = -1;

		List<Vertex> la = new LinkedList<Vertex>();
		la.add(v);

		clusters.put(all.getClusterSize(), la);

		// Add clusters of size 1.
		List<Vertex> l1 = new LinkedList<Vertex>();
		for (String t : stTaxa) {
			STITreeCluster c = new STITreeCluster(stTaxa);
			c.addLeaf(t);

			v = new Vertex();
			v._cluster = c;
			v._el_num = DeepCoalescencesCounter.getClusterCoalNum(trees, c, taxonMap, true);
			if(v._el_num > maxEL){
				maxEL = v._el_num;
			}
			v._min_cost = -1;

			l1.add(v);
		}
		maxEL = maxEL+1;
		clusters.put(1, l1);
		return maxEL;
	}

	/**
	 * This function is written to compute clusters induced by gene trees when there is single gene copy.
	 *
	 * @param 	trees	the list of gene trees
	 * @param	taxa	taxa in species tree/gene trees
	 * @param	clusters	clusters with the number of extra lineages
	 * @return  the maximal number of extra lineages
	 */
	private double computeTreeClusters(List<MutableTuple<Tree,Double>> trees, String taxa[], Map<Integer, List<Vertex>> clusters) {
		// Compute the list of vertices in the graph, which correspond to the induced clusters.
		double maxEL = 0;
        for(MutableTuple<Tree,Double> tr: trees){
			for (STITreeCluster tc : tr.Item1.getClusters(taxa, false)) {
				int tc_size = tc.getClusterSize();

				if (clusters.containsKey(tc_size)) {
					List<Vertex> l = clusters.get(tc_size);
					boolean found = false;

					for (Vertex v : l) {
						if (v._cluster.equals(tc)) {
							found = true;
							break;
						}
					}

					if (!found) {
						Vertex nv = new Vertex();
						nv._cluster = tc;
						nv._el_num = DeepCoalescencesCounter.getClusterCoalNum(trees, tc, true);
						if(nv._el_num > maxEL){
							maxEL = nv._el_num;
						}
						nv._min_cost = -1;

						l.add(nv);
					}
				}
				else {
					List<Vertex> l = new LinkedList<Vertex>();
					Vertex v = new Vertex();

					v._cluster = tc;
					v._el_num = DeepCoalescencesCounter.getClusterCoalNum(trees, tc, true);
					if(v._el_num > maxEL){
						maxEL = v._el_num;
					}
					v._min_cost = -1;

					l.add(v);
					clusters.put(tc_size, l);
				}
			}
		}

		// Add the cluster containing all taxa to the end of the list.
		STITreeCluster all = new STITreeCluster(taxa);
		for (String t : taxa) {
			all.addLeaf(t);
		}

		Vertex v = new Vertex();
		v._cluster = all;
		v._el_num = 0;
		v._min_cost = -1;

		List<Vertex> la = new LinkedList<Vertex>();
		la.add(v);

		clusters.put(all.getClusterSize(), la);

		// Add clusters of size 1.
		List<Vertex> l1 = new LinkedList<Vertex>();
		for (String t : taxa) {
			STITreeCluster c = new STITreeCluster(taxa);
			c.addLeaf(t);

			v = new Vertex();
			v._cluster = c;
			v._el_num = 0;
			v._min_cost = -1;

			l1.add(v);
		}
		maxEL = maxEL+1;
		clusters.put(1, l1);
		return maxEL;
	}

	/**
	 * Compute the minimum cost for trees whose leaf set is vertex's cluster.
	 *
	 * @param	clusters	clusters with the number of extra lineages
	 * @param	v	the vertex being computed
	 * @param	maxEL  the maximal number of extra lineages
	 * @return
	 */
	private double computeMinCost(Map<Integer, List<Vertex>> clusters, Vertex v, double maxEL) {
		if (v._max_score != -1) {
			return v._max_score;
		}

		if (v._cluster.getClusterSize() <= 1) {
			v._min_cost = v._el_num;
			v._max_score = maxEL - v._min_cost;
			v._min_lc = v._min_rc = null;

			return v._max_score;
		}

		int vsize = v._cluster.getClusterSize();

		for (int i = 1; i <= (vsize / 2); i++) {
			List<Vertex> leftList = clusters.get(i);
			if (leftList != null) {
				for (Vertex lv : leftList) {
					if (v._cluster.containsCluster(lv._cluster)) {
						double lscore = computeMinCost(clusters, lv, maxEL);

						List<Vertex> rightList = clusters.get(vsize - i);
						if (rightList != null) {
							for (Vertex rv : rightList) {
								if (lv._cluster.isDisjoint(rv._cluster) && v._cluster.containsCluster(rv._cluster)) {
									double rscore = computeMinCost(clusters, rv, maxEL);

									if ((v._max_score == -1) || (lscore + rscore + maxEL - v._el_num > v._max_score)) {
										v._min_cost = lv._min_cost + rv._min_cost + v._el_num;
										v._max_score = lscore + rscore + maxEL - v._el_num;
										v._min_lc = lv;
										v._min_rc = rv;
									}

									break;	// Already found the only pair of clusters whose union is v's cluster.
								}
							}
						}
					}
				}
			}
		}


		if (v._min_cost == -1) {	// Cannot resolve v as a binary node.
			STITreeCluster cl = v._cluster;
			ArrayList<Vertex> subClusters = new ArrayList<Vertex>();
			//get all the subclusters
			for(int i=1; i<cl.getClusterSize(); i++){
				List<Vertex> cls = clusters.get(i);
				if(cls == null)continue;
				for(Vertex tv:cls){
					if(cl.containsCluster(tv._cluster)){
						if(tv._max_score==-1){
							computeMinCost(clusters, tv, maxEL);
						}
						if(i>=2){
							subClusters.add(tv);
						}
					}
				}
			}

			for(int i=1;i<subClusters.size();i++){
				Vertex v1 = subClusters.get(i);
				for(int j=0;j<i;j++){
					Vertex v2 = subClusters.get(j);
					if(v1._el_num < v2._el_num){
						subClusters.remove(v1);
						subClusters.add(j,v1);
						break;
					}
					else if(v1._el_num == v2._el_num){
						if(v1._cluster.getClusterSize() > v2._cluster.getClusterSize()){
							subClusters.remove(v1);
							subClusters.add(j,v1);
							break;
						}
					}
				}
			}

/*
			for(int i=1;i<subClusters.size();i++){
				Vertex v1 = subClusters.get(i);
				for(int j=0;j<i;j++){
					Vertex v2 = subClusters.get(j);
					if(v1._cluster.getClusterSize() == v2._cluster.getClusterSize()){
						if(v1._max_score > v2._max_score){
							subClusters.remove(v1);
							subClusters.add(j,v1);
							break;
						}
					}
					else if(v1._cluster.getClusterSize() < v2._cluster.getClusterSize()){
						subClusters.remove(v1);
						subClusters.add(j,v1);
						break;
					}
				}
			}
*/


			ArrayList<Vertex> maxSubClusters = new ArrayList<Vertex>();
			List<Vertex> l1 = clusters.get(1);
			for (Vertex s : l1) {
				if(cl.containsCluster(s._cluster)){
					maxSubClusters.add(s);
				}
			}

			for(Vertex tv: subClusters){
				boolean compatible = true;
				boolean contained = false;
				ArrayList<Vertex> rmCls = new ArrayList<Vertex>();
				for(Vertex ex_v: maxSubClusters){
					if(ex_v._cluster.containsCluster(tv._cluster)){
						contained = true;
						break;
					}
					if(!tv._cluster.isCompatible(ex_v._cluster)){
						compatible = false;
						break;
					}
					else{
						if(tv._cluster.containsCluster(ex_v._cluster)){
							rmCls.add(ex_v);
						}
					}
				}
				if(compatible && (!contained)){
					maxSubClusters.removeAll(rmCls);
					maxSubClusters.add(tv);
				}
			}

			v._min_cost = v._el_num;
			v._max_score = maxEL - v._el_num;
			for(Vertex tv: maxSubClusters){
				v._min_cost += tv._min_cost;
				v._max_score += tv._max_score;
			}
			v._subcl = maxSubClusters;
		}

		return v._max_score;
	}



	/**
	 * Generate all possible clusters, not just those induced by gene trees.
	 *
	 * @param trees
	 * @param stTaxa
	 * @param taxonMap
	 * @param clusters
	 *
	 * @return the maximal number of clusters
	 */
	private double computeAllClusters(List<MutableTuple<Tree,Double>> trees, String stTaxa[], Map<String, String> taxonMap, Map<Integer, List<Vertex>> clusters) {
		int n = stTaxa.length;
		double maxEL = 0;
		if (n <= 0) {
			System.err.println("Empty list of taxa.");
			return -1;
		}

		BitSet counter = new BitSet(n);
		boolean done = false;

		while (!done) {	// Repeat until all 2^n - 1 binary strings are generated.
			int i = 0;
			while (i < n && counter.get(i)) {
				counter.clear(i);
				i++;
			}
			if (i >= n) {
				done = true;	// Already generated all binary strings.
			}
			else {
				counter.set(i, true);
				STITreeCluster tc = new STITreeCluster(stTaxa);
				tc.setCluster((BitSet) counter.clone());
				Vertex v = new Vertex();
				v._cluster = tc;
				if(taxonMap == null){
					v._el_num = DeepCoalescencesCounter.getClusterCoalNum(trees, tc, true);
				}
				else{
					v._el_num = DeepCoalescencesCounter.getClusterCoalNum(trees, tc, taxonMap, true);
				}
				if(v._el_num > maxEL){
					maxEL = v._el_num;
				}

				int size = tc.getClusterSize();
				List<Vertex> l = clusters.get(size);
				if (l == null) {
					l = new LinkedList<Vertex>();
				}

				l.add(v);
				clusters.put(size, l);
			}
		}
		maxEL = maxEL + 1;
		return maxEL;
	}

	/**
	 * Generate all possible clusters, not just those induced by gene trees.
	 *
	 * @param trees
	 * @param stTaxa
	 * @param clusters
	 *
	 * @return the maximal number of clusters
	 */
	private double computeAllClusters(List<MutableTuple<Tree,Double>> trees, String stTaxa[], Map<Integer, List<Vertex>> clusters) {
		return computeAllClusters(trees, stTaxa, null, clusters);
	}



	private int getResolutionsNumber(int nodeNumber){
		int total = 1;
		for(int i=3; i<=nodeNumber; i++){
			total = total * (2 * i - 3);
		}
		return total;
	}

	private double tryBinaryResolutions(Tree tr, double time, String[] taxa, List<MutableTuple<Tree,Double>> gts, Map<String,String> taxonMap){
		List<TNode> nodelist = new ArrayList<TNode>();
		List<Integer> degreelist = new ArrayList<Integer>();
		int totalResolutions = 0;
		for (TNode node : new PostTraversal<Object>(tr.getRoot())) {
			int childCount = node.getChildCount();
			if(childCount>2){
				nodelist.add(node);
				int resolutionsNumber = getResolutionsNumber(childCount);
				degreelist.add(resolutionsNumber);
				totalResolutions = totalResolutions + resolutionsNumber;
			}
		}
		double addedxl = 0;
		for(int i=0; i<nodelist.size(); i++){
			TNode unresolvedNode = nodelist.get(i);
			Map<Integer, TNode> id2node = new HashMap<Integer, TNode>();
			for(TNode child: unresolvedNode.getChildren()){
				id2node.put(child.getID(), child);
			}
			Integer[] childIDs = id2node.keySet().toArray(new Integer[0]);
			Map<BitSet, Double> cluster2xl = new HashMap<BitSet,Double>();
			double endtime;
			if(time == -1){
				endtime = -1;
			}
			else{
				endtime = time*(((double)degreelist.get(i))/totalResolutions)*1000 + System.currentTimeMillis();
			}
			Solution sol = addMoreLeaves(null , childIDs, 0, id2node, taxa, gts, taxonMap, endtime , cluster2xl);
			TNode parent = unresolvedNode.getParent();
			double xl = ((STINode<Double>)unresolvedNode).getData();
			((STINode)unresolvedNode).removeNode();
			if(parent!=null){
				TNode newnode = ((STINode)parent).createChild(sol.getTree().getRoot());
				((STINode<Double>)newnode).setData(xl);
			}
			else{
				for(TNode child: sol.getTree().getRoot().getChildren()){
					((STINode)tr.getRoot()).createChild(child);
				}
			}
			addedxl = addedxl + sol.getCoalNum();
		}
		return addedxl;
	}


	private Solution addMoreLeaves(STITree<Integer> preTree, Integer[] leavesid, int index, Map<Integer, TNode> id2node, String[] taxa, List<MutableTuple<Tree,Double>> gts, Map<String,String> taxonMap, double endTime, Map<BitSet, Double> cluster2xl){
		if(preTree == null){
			preTree = new STITree<Integer>(false);
			STINode<Integer> root = preTree.getRoot();
			STINode<Integer> newnode = root.createChild();
			newnode.setData(leavesid[index++]);
			STINode<Integer> innode = root.createChild();
			newnode = innode.createChild();
			newnode.setData(leavesid[index++]);
			newnode = innode.createChild();
			newnode.setData(leavesid[index++]);
		}
		Solution sol = null;
		if(index == leavesid.length){
			sol = tryAllRootings(preTree, id2node, taxa, gts, taxonMap, endTime, cluster2xl);
		}
		else{
			int id = leavesid[index++];
			for(TNode n: new PostTraversal<Object>(preTree.getRoot())){
				if(!n.isLeaf()){
					if(n.getChildCount()!=2){
						throw new RuntimeException("Not binary!");
					}
					Iterator it = n.getChildren().iterator();
					for(int i=0; i<2; i++){
						STINode<Integer> child = (STINode<Integer>)(it.next());
						STITree<Integer> newTree = new STITree<Integer>(preTree);
						TNode peerChild = newTree.getNode(child.getID());
						TNode peerParent = peerChild.getParent();
						STINode<Integer> newchild = ((STINode<Integer>)peerParent).createChild();
						newchild.adoptChild((TMutableNode)peerChild);
						STINode<Integer> newnode = newchild.createChild();
						newnode.setData(id);
						Solution thissol = addMoreLeaves(newTree, leavesid, index, id2node, taxa, gts, taxonMap, endTime, cluster2xl);
						if(sol==null || sol.getCoalNum()>thissol.getCoalNum()){
							sol = thissol;
						}
						double now = System.currentTimeMillis();
						if(endTime != -1 && now > endTime){
							return sol;
						}
						if(n.isRoot()){
							break;
						}
					}
				}
			}
		}
		return sol;
	}

	private Solution tryAllRootings(Tree subtree, Map<Integer, TNode> id2node, String[] taxa, List<MutableTuple<Tree,Double>> gts, Map<String,String> taxonMap, double endTime, Map<BitSet, Double> cluster2xl){
		Solution sol = new Solution();
		sol._totalCoals = -1;
		for(Tree rootedsubtree: subtree.getAllRootingTrees()){
			for(TNode replacingNode : new PostTraversal<Object>(rootedsubtree.getRoot())){
				if(replacingNode.isLeaf()){
					int id = ((STINode<Integer>)replacingNode).getData();
					TNode peerNode = id2node.get(id);
					if(peerNode.isLeaf()){
						((STINode)replacingNode).setName(peerNode.getName());
						//STINode<Integer> newNode = ((STINode)replacingNode).createChild(peerNode.getName());
						//newNode.setData(((STINode<Integer>)peerNode).getData());
					}
					else{
						for(TNode child: peerNode.getChildren()){
							((STINode)replacingNode).createChild(child);
						}
						((STINode)replacingNode).setName("");
						//((STINode<Integer>)replacingNode).setData(((STINode<Integer>)peerNode).getData());
					}
					((STINode<Integer>)replacingNode).setData(((STINode<Integer>)peerNode).getData());
				}
			}
			//Trees.removeBinaryNodes((MutableTree)rootedsubtree);
			//System.out.println(rootedsubtree.toStringWD());
			double xl = 0;
			Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
			for (TNode node : new PostTraversal<Integer>(rootedsubtree.getRoot())) {
				BitSet bs = new BitSet();
				if (node.isLeaf()) {
					for (int i = 0; i < taxa.length; i++) {
						if (node.getName().equals(taxa[i])) {
							bs.set(i);
							break;
						}
					}
					map.put(node, bs);
				}
				else {
					for (TNode child : node.getChildren()) {
						BitSet childCluster = map.get(child);
						bs.or(childCluster);
					}
					map.put(node, bs);
				}
				if(!node.isRoot()){
					Double el = ((STINode<Double>)node).getData();
					if(el == null){
						el = cluster2xl.get(bs);
						if(el == null){
							STITreeCluster tc = new STITreeCluster(taxa);
							tc.setCluster(bs);
							if(taxonMap==null){
								el = DeepCoalescencesCounter.getClusterCoalNum(gts, tc, true);
							}
							else{
								el = DeepCoalescencesCounter.getClusterCoalNum(gts, tc, taxonMap, true);
							}
							cluster2xl.put(bs, el);
						}
						((STINode<Double>)node).setData(el);
						xl = xl + el;
					}
				}
			}
			//System.out.println(rootedSubtree);
			if(sol.getCoalNum()==-1 || sol.getCoalNum()>xl){
				sol._totalCoals = xl;
				sol._st = rootedsubtree;
			}
		}
		return sol;
	}

	class Vertex {
		public STITreeCluster _cluster;	// The cluster associated with this class.
		public double _el_num;				// Number of extra lineages for this cluster.
		public double _min_cost;			// Minimum cost in all trees whose leaf set is _cluster.
		public double _max_score;			// The number of score being used for finding the optimal tree in dynamic programming
		public Vertex _min_lc;			// Left cluster for the min tree.
		public Vertex _min_rc;			// Right cluster for the min tree.
		public List<Vertex> _subcl;		// The list of child clusters which is used when this cannot be solved as a binary node

		public Vertex() {
			_cluster = null;
			_el_num = _min_cost = _max_score = -1;
			_min_lc = _min_rc = null;
			_subcl = null;
		}

		public String toString(){
			return _cluster.toString()+"/"+_el_num;
		}
	}
}
