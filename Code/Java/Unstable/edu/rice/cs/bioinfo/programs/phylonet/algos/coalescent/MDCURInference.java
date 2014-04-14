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

import edu.rice.cs.bioinfo.programs.phylonet.algos.MaxClique;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/10/11
 * Time: 1:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class MDCURInference
{


	/**
	 * Return the minimum number of extra lineages for the optimal tree. This function should be called
	 * only after the optimal tree is already inferred (by calling inferSpeciesTree).
	 *
	 * @return
	 */
/*	public int getMinExtraLineageNumber() {
		if (_min_el_num >= 0) {
			return _min_el_num;
		}
		else {
			System.err.println("The optimal tree has not been inferred. The function returns -1.");
			return -1;
		}
	}
*/
	/**
	 * Return the total number of clusters induced by all gene trees. This function should be called only after
	 * the optimal tree is already inferred (by calling inferSpeciesTree).
	 *
	 * @return
	 */
/*
	public int getTotalClusterNum() {
		if (clusters == null) {
			System.err.println("The optimal tree has not been inferred. The function returns -1.");
			return -1;
		}
		else {
			int total = 0;
			for (List<Vertex> l : clusters.values()) {
				total += l.size();
			}

			return total - 1;	// Excluding the all-taxon cluster.
		}
	}
*/
	/**
	 * Infers the species tree from the given list of unrooted gene trees with single allele.
	 *
	 * @param 	trees	the list of gene trees from which the species tree is to be inferred
	 * @param	explore	false if the result is the only one optimal tree
	 *                  true if the result is a set of trees
	 * @param   proportion	  specified only when explore is set to true
	 *          			  specifies the range of the number of extra lineages that the set of trees have
	 * @param   exhaust	  false if clusters are induced from gene trees
	 * 					  true if clusters all all possible clusters
	 * @return	tree(s) inferred from the supplied gene trees with the total number of extra lineages
	 */
	public List<Solution> inferSpeciesTree(List<Tree> trees, boolean explore, double proportion, boolean exhaust) {
		if (trees == null || trees.size() == 0) {
			throw new IllegalArgumentException("empty or null list of trees");
		}

		String taxa[] = new String[trees.get(0).getLeafCount()];
		int count = 0;
		for (TNode node : trees.get(0).postTraverse()) {
			if (node.isLeaf()) {
				taxa[count++] = node.getName();
			}
		}

		// Find the tree with minimum score over all taxa.
		Map<Integer, List<Vertex>> clusters;
		if(!exhaust){
			clusters = computeTreeClusters(trees, taxa);
		}
		else{
			clusters = computeAllClusters(trees, taxa);
		}
		List<Solution> solutions;
		if(explore){
			solutions = findTreesByClique(clusters,taxa,proportion);
		}
		else{
			solutions = findTreesByDP(clusters,taxa);
		}
		return solutions;

	}

	/**
	 * Infers the species tree from the given list of unrooted gene trees with multiple alleles.
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
	public List<Solution> inferSpeciesTree(List<Tree> trees, Map<String, String> taxonMap, boolean explore, double proportion, boolean exhaust) {
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
		Map<Integer, List<Vertex>> clusters;
		if(!exhaust){
			clusters = computeTreeClusters(trees, stTaxa, gtTaxa, taxonMap);
		}
		else{
			clusters = computeAllClusters(trees, stTaxa, gtTaxa, taxonMap);
		}
		List<Solution> solutions;
		if(explore){
			solutions = findTreesByClique(clusters,stTaxa,proportion);
		}
		else{
			solutions = findTreesByDP(clusters,stTaxa);
		}
		return solutions;
	}

	/**
	 * Find the species tree by dynamic programming
	 *
	 * @param 	clusters	clusters with extra lineages
	 * @param	stTaxa	  taxa expression
	 *
	 * @return	the tree inferred from the supplied gene trees with the total number of extra lineages
	 *
	 */
	private List<Solution> findTreesByDP(Map<Integer, List<Vertex>> clusters, String[] stTaxa){
		List<Solution> solutions = new ArrayList<Solution>();

		int maxEL = 0;
		for(Integer size: clusters.keySet()){
			if(size > stTaxa.length/2){
				//System.out.println(clusters.get(size));
				continue;
			}
			for(Vertex v:clusters.get(size)){
				if(v._c_el_num!=-1){
					continue;
				}
				STITreeCluster cc = v._cluster.complementaryCluster();
				for(Vertex cv: clusters.get(stTaxa.length-size)){
					if(cc.equals(cv._cluster)){
						if(v._c_el_num > maxEL){
							maxEL = v._c_el_num;
						}
						if(cv._el_num > maxEL){
							maxEL = cv._el_num;
						}
						v._c_el_num = cv._el_num;
						cv._c_el_num = v._el_num;
						break;
					}
				}
			}
			//System.out.println(clusters.get(size));
		}

		maxEL ++;

		Vertex all = clusters.get(stTaxa.length).get(0);
		computeMinCost(clusters,all,maxEL);
		// Build the minimum tree.
		List<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
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
	 *
	 * @return	the tree inferred from the supplied gene trees with the total number of extra lineages
	 *
	 */
	private List<Solution> findTreesByClique(Map<Integer, List<Vertex>> cmap, String[] stTaxa, double proportion){
		List<Solution> solutions = new LinkedList<Solution>();

		List<STITreeCluster<Integer>> clusters = new LinkedList<STITreeCluster<Integer>>();
		int addEL = 0;

		for(Map.Entry<Integer, List<Vertex>> entry: cmap.entrySet()){
			if(entry.getKey()==1 || entry.getKey()==stTaxa.length-1){
				List<Vertex> l = entry.getValue();
				for(Vertex v: l){
					addEL += v._el_num;
				}
			}
			if(entry.getKey() < stTaxa.length && entry.getKey()>1){
				List<Vertex> l = entry.getValue();
				for(Vertex v: l){
					STITreeCluster<Integer> c = new STITreeCluster<Integer>(v._cluster);
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
				List<Integer> cls = new ArrayList<Integer>();
				for(int id: nodes){
					STITreeCluster<Integer> c = clusters.get(id);
					if(c.getClusterSize()!=stTaxa.length-1 && !cls.contains(id)){

						sum += c.getData();
						STITreeCluster cc = c.complementaryCluster();
						int ccID= clusters.indexOf(cc);
						sum += clusters.get(ccID).getData();
						cls.add(id);
						cls.add(ccID);
					}
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
		SymmetricDifference sd = new SymmetricDifference();
		for(Solution s: maxCliques){
			if(s._totalCoals <= maxCoal){
				List<STITreeCluster> minClusters = new ArrayList<STITreeCluster>();
				for(int id: s._clusterIDs){
					minClusters.add(clusters.get(id));
				}
				s._st = Trees.buildTreeFromClusters(minClusters);
				boolean dup = false;
				for(Solution ex_s: solutions){
					sd.computeDifference(s._st, ex_s._st, true);
					if(sd.getWeightedAverage()==0){
						dup = true;
						break;
					}
				}
				if(!dup){
					solutions.add(s);
				}
			}
			else{
				break;
			}
		}
		return solutions;
	}



	/**
	 * This function is written to compute clusters induced by gene trees when there are multiple
	 * gene copies.
	 *
	 * @param 	trees	the list of gene trees
	 * @param	stTaxa	taxa in species tree
	 * @param	gtTaxa	taxa in gene trees
	 * @param	taxonMap	Maps gene tree taxa to species tree taxa.
	 *
	 * @return  clusters with the number of extra lineages
	 */
	private Map<Integer, List<Vertex>> computeTreeClusters(List<Tree> trees, String stTaxa[], String gtTaxa[], Map<String, String> taxonMap) {
		Map<Integer, List<Vertex>> clusters = new HashMap<Integer, List<Vertex>>();

		List<List<STITreeCluster>> treeCls = new ArrayList<List<STITreeCluster>>();
		List<STITreeCluster> allCls = new LinkedList<STITreeCluster>();

		// Compute all the clusters

		for (Tree tr : trees) {
			List<STITreeCluster> treeCl = new ArrayList<STITreeCluster>();
			for(STITreeCluster gttc: tr.getBipartitionClusters(gtTaxa,false)){
				STITreeCluster tc = new STITreeCluster(stTaxa);
				for(String leaf: gttc.getClusterLeaves()){
					tc.addLeaf(taxonMap.get(leaf));
				}
				int size= gttc.getClusterSize();
				if(size<gtTaxa.length && size>1){
					if(!treeCl.contains(gttc)){
						treeCl.add(gttc);
					}
				}

				size= tc.getClusterSize();
				if(size<stTaxa.length){
					if(!allCls.contains(tc)){
						allCls.add(tc);
					}
				}
			}
			treeCls.add(treeCl);
		}

		int clsize = allCls.size();
		for(int i=0; i<clsize; i++){
			STITreeCluster cl = allCls.get(i);
			STITreeCluster clc = cl.complementaryCluster();
			if(!allCls.contains(clc)){
				allCls.add(clc);
			}
		}

		//compute all the vertex in graph
		for(STITreeCluster tc:allCls){
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
					nv._el_num = getClusterCoalNum(trees, tc, treeCls, gtTaxa, taxonMap);
					nv._min_cost = -1;

					l.add(nv);
				}
			}
			else {
				List<Vertex> l = new LinkedList<Vertex>();

				Vertex v = new Vertex();
				v._cluster = tc;
				v._el_num = getClusterCoalNum(trees, tc, treeCls, gtTaxa, taxonMap);
				v._min_cost = -1;

				l.add(v);
				clusters.put(tc_size, l);
			}
		}

		// Add the cluster containing all taxa to the end of the list.
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

		return clusters;
	}

	/**
	 * This function is written to compute clusters induced by gene trees when there are multiple
	 * gene copies.
	 *
	 * @param 	trees	the list of gene trees
	 * @param	taxa	taxa in species tree/gene trees
	 *
	 * @return  clusters with the number of extra lineages
	 */
	private Map<Integer, List<Vertex>> computeTreeClusters(List<Tree> trees, String taxa[]) {
		Map<Integer, List<Vertex>> clusters = new HashMap<Integer, List<Vertex>>();

		List<List<STITreeCluster>> treeCls = new ArrayList<List<STITreeCluster>>();
		List<STITreeCluster> allCls = new LinkedList<STITreeCluster>();

		// Compute all the clusters

		for (Tree tr : trees) {
			List<STITreeCluster> treeCl = new ArrayList<STITreeCluster>();
			for(STITreeCluster tc:tr.getBipartitionClusters(taxa,false)){
				int tc_size= tc.getClusterSize();
				if(tc_size<taxa.length){
					if(tc_size > 1){
						if(!treeCl.contains(tc)){
							treeCl.add(tc);
						}
					}
					if(!allCls.contains(tc)){
						allCls.add(tc);
					}
				}
			}
			treeCls.add(treeCl);
		}

		//compute all the vertex in graph
		for(STITreeCluster tc:allCls){
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
					if(tc.getClusterSize()==1){
						nv._el_num = 0;
					}
					else{
						nv._el_num = getClusterCoalNum(trees, tc, treeCls);
					}
					nv._min_cost = -1;

					l.add(nv);
				}
			}
			else {
				List<Vertex> l = new LinkedList<Vertex>();

				Vertex v = new Vertex();
				v._cluster = tc;
				if(tc.getClusterSize()==1){
					v._el_num = 0;
				}
				else{
					v._el_num = getClusterCoalNum(trees, tc, treeCls);
				}
				v._min_cost = -1;
				l.add(v);
				clusters.put(tc_size, l);
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

		return clusters;
	}

	/**
	 * Generate all possible clusters, not just those induced by gene trees.
	 *
	 * @param trees
	 * @param stTaxa
	 *
	 * @return resulting clusters
	 */
	private Map<Integer, List<Vertex>> computeAllClusters(List<Tree> trees, String stTaxa[]) {
		return computeAllClusters(trees, stTaxa, null, null);
	}

	/**
	 * Generate all possible clusters, not just those induced by gene trees.
	 *
	 * @param trees
	 * @param stTaxa
	 * @param gtTaxa
	 * @param taxonMap
	 *
	 * @return resulting clusters
	 */
	private Map<Integer, List<Vertex>> computeAllClusters(List<Tree> trees, String stTaxa[], String gtTaxa[], Map<String, String> taxonMap) {
		int n = stTaxa.length;
		if (n <= 0) {
			System.err.println("Empty list of taxa.");
			return null;
		}
		Map<Integer, List<Vertex>> clusters = new HashMap<Integer, List<Vertex>>();
		// Compute all the clusters in trees
		List<List<STITreeCluster>> treeCls = new ArrayList<List<STITreeCluster>>();
		for (Tree tr : trees) {
			String[] taxa;
			if(taxonMap == null){
				taxa = stTaxa;
			}
			else{
				taxa = gtTaxa;
			}
			List<STITreeCluster> treeCl = new ArrayList<STITreeCluster>();
			for(STITreeCluster tc: tr.getBipartitionClusters(taxa,false)){
				int size= tc.getClusterSize();
				if(size<taxa.length && size>1){
					if(!treeCl.contains(tc)){
						treeCl.add(tc);
					}
				}
			}
			treeCls.add(treeCl);
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
				if(tc.getClusterSize() == stTaxa.length){
					v._el_num = 0;
					v._c_el_num = 0;
				}
				else {
					if(taxonMap == null){
						v._el_num = getClusterCoalNum(trees, tc, treeCls);
					}
					else{
						v._el_num = getClusterCoalNum(trees, tc, treeCls, gtTaxa, taxonMap);
					}
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
		return clusters;
	}



	/**
	 * Compute the minimum cost for trees whose leaf set is vertex's cluster.
	 *
	 * @param	clusters	clusters with the number of extra lineages
	 * @param	v	the vertex being computed
	 * @param	maxEL  the maximal number of extra lineages
	 * @return
	 */
	private int computeMinCost(Map<Integer, List<Vertex>> clusters, Vertex v, int maxEL) {
		if (v._max_score != -1) {
			return v._max_score;
		}

		if (v._cluster.getClusterSize() <= 1) {
			v._min_cost = v._el_num + v._c_el_num;
			v._max_score = 2*maxEL - v._min_cost;
			v._min_lc = v._min_rc = null;
			return v._max_score;
		}

		int vsize = v._cluster.getClusterSize();

		for (int i = 1; i <= (vsize / 2); i++) {
			List<Vertex> leftList = clusters.get(i);
			if (leftList != null) {
				for (Vertex lv : leftList) {
					if (v._cluster.containsCluster(lv._cluster)) {
						int lscore = computeMinCost(clusters, lv, maxEL);
						List<Vertex> rightList = clusters.get(vsize - i);
						if (rightList != null) {
							for (Vertex rv : rightList) {
								if (lv._cluster.isDisjoint(rv._cluster) && v._cluster.containsCluster(rv._cluster)) {
									int rscore = computeMinCost(clusters, rv, maxEL);
									int score = 0 ,cost = 0;
									if(vsize!=v._cluster.getTaxa().length){
										score = lscore + rscore + 2*maxEL - v._el_num - v._c_el_num;
										cost = lv._min_cost + rv._min_cost + v._el_num + v._c_el_num;
									}
									else{
										score = lscore + rscore + maxEL - v._el_num - (2*maxEL - (rv._c_el_num + rv._el_num));
										cost = lv._min_cost + rv._min_cost + v._el_num - rv._c_el_num - rv._el_num;
									}
									if ((v._max_score == -1) || (score > v._max_score)) {
										v._min_cost = cost;
										v._max_score = score;
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
			//System.out.println("yes");
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
					for(Vertex ex_v: rmCls){
						maxSubClusters.remove(ex_v);
					}
					maxSubClusters.add(tv);
				}
			}
			//System.out.println(v);
			v._min_cost = v._el_num + v._c_el_num;
			v._max_score = 2*maxEL - v._min_cost;
			for(Vertex tv:maxSubClusters){
			//	System.out.println(tv+" "+tv._min_cost);
				v._min_cost = v._min_cost + tv._min_cost;
				v._max_score = v._max_score + tv._max_score;
			}
			//System.out.println();
			v._subcl = maxSubClusters;

		}


	/*

		if (v._min_cost == -1) {	// Cannot resolve v as a binary node.
			// Consider this as a big star.
			v._min_cost = v._el_num + v._c_el_num;
			v._max_score = 2*maxEL - v._min_cost;
			List<Vertex> l1 = clusters.get(1);
			for (Vertex s : l1) {
				if (v._cluster.containsCluster(s._cluster)) {
					if(s._max_score == -1){
						computeMinCost(clusters,s,maxEL);
					}
					v._min_cost += s._min_cost;
					v._max_score += s._max_score;

				}
			}
		}
*/

		return v._max_score;
	}


	/**
	 * Compute the number of extra lineages contributed by a cluster in a set of trees when single allele is involved
	 *
	 * @param trees		a list of given gene trees
	 * @param cluster	the cluster being computed
	 * @param treeCls	a list of a list of clusters induced by each gene tree
	 *
	 * @return the number of extra lineage
	 */
	public int getClusterCoalNum(List<Tree> trees, STITreeCluster cluster, List<List<STITreeCluster>> treeCls) {
		int weight = 0;

		int index = 0;
		for (Tree tr : trees) {

			weight += getClusterCoalNum(cluster, treeCls.get(index++));
		}

		return weight;
	}

	/**
	 * Compute the number of extra lineages contributed by a cluster in a set of trees when multiple alleles are involved
	 *
	 * @param trees		a list of given gene trees
	 * @param cluster	the cluster being computed
	 * @param treeCls	a list of a list of clusters induced by each gene tree
	 * @param gtTaxa	taxa of gene trees
	 * @param taxonMap	Maps gene tree taxa to species tree taxa.
	 *
	 * @return the number of extra lineage
	 */
	public int getClusterCoalNum(List<Tree> trees, STITreeCluster cluster, List<List<STITreeCluster>> treeCls, String[] gtTaxa,Map<String,String> taxonMap) {
		int weight = 0;

		int index = 0;
		for (Tree tr : trees) {
			STITreeCluster gtcl = new STITreeCluster(gtTaxa);
			for(TNode n: tr.getNodes()){
				if(n.isLeaf()){
					if(cluster.containsLeaf(taxonMap.get(n.getName()))){
						gtcl.addLeaf(n.getName());
					}
				}
			}
			weight += getClusterCoalNum(gtcl, treeCls.get(index++));
		}

		return weight;
	}



	/**
	 * Compute the number of extra lineages contributed by a cluster in one tree
	 *
	 * @param cluster	the cluster being computed
	 * @param treeCl	a list of clusters induced by the given gene tree
	 *
	 * @return the number of extra lineage
	 */
	public int getClusterCoalNum(STITreeCluster cluster, List<STITreeCluster> treeCl) {
		List<STITreeCluster> maxSubClusters = new ArrayList<STITreeCluster>();
		String[] taxa = cluster.getTaxa();
		for(String leaf: cluster.getClusterLeaves()){
			STITreeCluster newCl = new STITreeCluster(taxa);
			newCl.addLeaf(leaf);
			maxSubClusters.add(newCl);
		}


		for(STITreeCluster cl: treeCl){
			if(cluster.containsCluster(cl)){
				List<STITreeCluster> rmCls = new ArrayList<STITreeCluster>();
				for(STITreeCluster maxSubCl: maxSubClusters){
					if(cl.containsCluster(maxSubCl)){
						rmCls.add(maxSubCl);
					}
				}
				if(rmCls.size()>0){
					for(STITreeCluster rmcl: rmCls){
						maxSubClusters.remove(rmcl);
					}
					maxSubClusters.add(cl);
				}
				if(maxSubClusters.size()==1){
					break;
				}
			}
		}
		return maxSubClusters.size()-1;
	}
	class Vertex {
		public STITreeCluster _cluster;	// The cluster associated with this class.
		public int _el_num;				// Number of extra lineages for this cluster.
		public int _c_el_num;			// Number of extra lineages for the complementary cluster
		public int _min_cost;			// Minimum cost in all trees whose leaf set is _cluster.
		public int _max_score;			// The number of score being used for finding the optimal tree in dynamic programming
		public Vertex _min_lc;			// Left cluster for the min tree.
		public Vertex _min_rc;			// Right cluster for the min tree.
		public List<Vertex> _subcl;		// The list of child clusters which is used when this cannot be solved as a binary node

		public Vertex() {
			_cluster = null;
			_el_num = _c_el_num= _min_cost = _max_score= -1;
			_min_lc = _min_rc = null;
			_subcl = null;
		}


		public String toString(){
			return _cluster.toString()+"/"+_el_num+"/"+_c_el_num;
		}
	}
}
