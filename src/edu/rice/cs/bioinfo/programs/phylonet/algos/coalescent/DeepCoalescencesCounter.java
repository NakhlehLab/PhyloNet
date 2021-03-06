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
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by Yun Yu
 * Date: 10/17/11
 * Time: 3:24 PM
 *
 * This class is to compute the minimum number of extra lineages required to reconcile gene trees within a species tree
 * For all functions related to rooted gene trees, see "Species tree inference by minimizing deep coalescences", PLoS Computational Biology, 2009 for details.
 * For all functions related to unrooted gene trees, see "Algorithms for MDC-based multi-locus phylogeny inference: Beyond rooted binary gene trees on single alleles", Journal of Computational Biology, 2011 for details.
 */
public class DeepCoalescencesCounter {

	/**
	 * Computes the minimum number of extra lineages given a species tree and a set of gene trees when single allele is sampled per species
	 *
	 * @param gts	a list of given gene trees
	 * @param st	a given species tree
	 * @param rooted	whether treating the gene trees as rooted
	 * @param bootstrap	threshold of bootstrap value. Branches with lower bootstrap values will be contracted
	 *
	 * @return the number of extra lineage
	 */
	public static double countExtraCoal(List<MutableTuple<Tree,Double>> gts,Tree st, boolean rooted, double bootstrap){
		double sum = 0;
		String[] taxa = st.getLeaves();

		if(bootstrap<100){
			for(MutableTuple<Tree,Double> tr: gts){
				if(Trees.handleBootStrapInTree(tr.Item1, bootstrap)==-1){
					throw new IllegalArgumentException("Input gene trees have nodes that don't have bootstrap value");
				}
			}
		}

		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		for (TNode node : st.postTraverse()) {
			BitSet bs = new BitSet();
			if (node.isLeaf()) {
				for (int i = 0; i < taxa.length; i++) {
					if (node.getName().equals(taxa[i])) {
						bs.set(i);
						break;
					}
				}
				map.put(node, bs);
				((STINode<Double>)node).setData(0.0);
			}
			else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = map.get(child);
					bs.or(childCluster);
				}
				map.put(node, bs);
				STITreeCluster c = new STITreeCluster(taxa);
				c.setCluster(bs);
				if(c.getClusterSize()==taxa.length){
					((STINode<Double>)node).setData(0.0);
				}
				else{
					double el = getClusterCoalNum(gts, c, rooted);
					((STINode<Double>)node).setData(el);
					sum += el;
				}
			}
		}
		return sum;
	}


	/**
	 * Computes the minimum number of extra lineages given a species tree and a set of gene tree when multiple alleles are sampled per species
	 *
	 * @param gts	a list of given gene trees
	 * @param st	a given species tree
	 * @param taxonMap	Maps gene tree taxa to species tree taxa
	 * @param rooted	whether treating the gene trees as rooted
	 * @param bootstrap	threshold of bootstrap value. Branches with lower bootstrap values will be contracted
	 *
	 * @return the number of extra lineage
	 */
	public static double countExtraCoal(List<MutableTuple<Tree,Double>> gts,Tree st, Map<String, String> taxonMap, boolean rooted, double bootstrap){
		String error = Trees.checkMapping(gts, taxonMap);
		if(error!=null){
			throw new RuntimeException("Gene trees have leaf named " + error + " that hasn't been defined in the mapping file");
		}

		int sum = 0;
		String[] stTaxa = st.getLeaves();

			if(bootstrap<1){
				for(MutableTuple<Tree,Double> tr: gts){
					if(Trees.handleBootStrapInTree(tr.Item1, bootstrap)==-1){
						throw new IllegalArgumentException("Input gene trees have nodes that don't have bootstrap value");
					}
				}
			}


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
				if(c.getClusterSize()==stTaxa.length){
					((STINode<Double>)node).setData(0.0);
				}
				else{
					double el = getClusterCoalNum(gts, c, taxonMap, rooted);
					((STINode<Double>)node).setData(el);
					sum += el;
				}
			}


		return sum;
	}



	/**
	 * Computes the minimum number of extra lineages given a cluster (in species tree) and a set of gene trees when single allele is sampled per species
	 *
	 * @param trees		a list of given gene trees
	 * @param cluster	the cluster being computed
	 * @param rooted	whether treating the gene trees as rooted
	 *
	 * @return the number of extra lineage
	 */
	public static double getClusterCoalNum(List<MutableTuple<Tree,Double>> trees, STITreeCluster cluster, boolean rooted) {
		double xl = 0;

        for(MutableTuple<Tree,Double> tr: trees){
			if(rooted){
				xl += getClusterCoalNum_rooted(tr.Item1, cluster) * tr.Item2;

			}
			else{
				xl += getClusterCoalNum_unrooted(tr.Item1, cluster) * tr.Item2;
			}
		}

		return xl;
	}

	/**
	 * Computes the minimum number of extra lineages given a cluster (in species tree) and a set of gene trees when multiple alleles are sampled per species
	 *
	 * @param trees		a list of given gene trees
	 * @param cluster	the cluster being computed
	 * @param taxonMap	Maps gene tree taxa to species tree taxa.
	 * @param rooted	whether treating the gene trees as rooted
	 *
	 * @return the number of extra lineage
	 */
	public static double getClusterCoalNum(List<MutableTuple<Tree,Double>> trees, STITreeCluster cluster, Map<String, String> taxonMap, boolean rooted) {
		double xl = 0;

		for (MutableTuple<Tree,Double> tr : trees) {
			if(rooted){
                xl += getClusterCoalNum_rooted(tr.Item1, cluster, taxonMap) * tr.Item2;
			}
			else{
                xl += getClusterCoalNum_unrooted(tr.Item1, cluster, taxonMap) * tr.Item2;
			}
		}

		return xl;
	}


	/**
	 * Computes the number of extra lineages given a cluster (in a species tree) and a rooted gene tree when single allele is sampled per species
	 *
	 * @param tr		a rooted gene tree
	 * @param cluster	the cluster being computed
	 *
	 * @return the number of extra lineage
	 */
	public static int getClusterCoalNum_rooted(Tree tr, STITreeCluster cluster) {
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		List<String> taxa = new LinkedList<String>();

		for (String t : cluster.getTaxa()) {
			taxa.add(t);
		}

		int count = 0;
		for (TNode node : tr.postTraverse()) {
			if (node.isLeaf()) {
				int index = taxa.indexOf(node.getName());
				BitSet bs = new BitSet();

				bs.set(index);
				if (cluster.containsCluster(bs)) {
					count++;
				}

				map.put(node, bs);
			}
			else {
				BitSet bs = new BitSet();
				int intersect = 0;
				int childCount = node.getChildCount();
				for (TNode child : node.getChildren()) {
					BitSet v = map.get(child);
					bs.or(v);
					if(childCount>2){
						if(cluster.containsCluster(v)){
							intersect ++;
						}
					}
				}

				if (cluster.containsCluster(bs)) {
					count -= node.getChildCount();
					count++;
				}
				else if(intersect>1){
					count -= intersect;
					count ++;
				}

				map.put(node, bs);
			}
		}
		return Math.max(count-1, 0);
	}


	/**
	 * Computes the number of extra lineages given a cluster (in a species tree) and an unrooted gene tree when single allele is sampled per species
	 *
	 * @param tr		a unrooted gene tree
	 * @param cluster	the cluster being computed
	 *
	 * @return the number of extra lineage
	 */
	public static int getClusterCoalNum_unrooted(Tree tr, STITreeCluster cluster) {
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		List<String> taxalist = new ArrayList<String>();
		String[] taxa = tr.getLeaves();
		int ntaxa = taxa.length;
		for(String leaf: taxa){
			taxalist.add(leaf);
		}
		STITreeCluster concluster = new STITreeCluster(taxa);
		for(String leaf:cluster.getClusterLeaves()){
			if(taxalist.contains(leaf)){
				concluster.addLeaf(leaf);
			}
		}
		if(concluster.getClusterSize()==ntaxa){
			return 0;
		}
		List<BitSet> coveragelist = new ArrayList<BitSet>();
		for (int i = concluster.getCluster().nextSetBit(0); i >= 0; i = concluster.getCluster().nextSetBit(i+1)) {
			BitSet bs = new BitSet(ntaxa);
			bs.set(i);
			coveragelist.add(bs);
		}
		for (TNode node :tr.postTraverse()) {
			if(coveragelist.size()<=1){
				break;
			}
			BitSet bs = new BitSet(ntaxa);
			int intersect = 0;
			BitSet virtualbs = new BitSet(ntaxa);
			if (node.isLeaf()) {
				int index = taxalist.indexOf(node.getName());
				bs.set(index);
				map.put(node, bs);
			}
			else {
				for (TNode child : node.getChildren()) {
					BitSet v = map.get(child);
					bs.or(v);
					if(concluster.containsCluster(v)){
						intersect ++;
						virtualbs.or(v);
					}
				}

				if (concluster.containsCluster(bs) || intersect>1) {
					if(concluster.containsCluster(bs)){
						virtualbs = bs;
					}
					for(int i=0; i<coveragelist.size();i++){
						BitSet exbs = coveragelist.get(i);
						BitSet temp = (BitSet) virtualbs.clone();
						temp.and(exbs);
						if(temp.equals(exbs)){
							coveragelist.remove(i);
							i--;
						}
					}
					coveragelist.add(virtualbs);
				}

				map.put(node, bs);
			}

			if(!node.isRoot()){
				BitSet complementbs = (BitSet)bs.clone();
				complementbs.flip(0, taxa.length);
				if(intersect>0){
					complementbs.or(virtualbs);
				}
				if(concluster.containsCluster(complementbs)){
					for(int i=0; i<coveragelist.size();i++){
						BitSet exbs = coveragelist.get(i);
						BitSet temp = (BitSet) complementbs.clone();
						temp.and(exbs);
						if(temp.equals(exbs)){
							coveragelist.remove(i);
							i--;
						}
					}
					coveragelist.add(complementbs);
					break;
				}
			}
		}
		/*
		if(coveragelist.size()!=1){
		System.out.println(cluster + ":" + tr+ ":"+ (coveragelist.size()-1));
		}
		*/
		return Math.max(0,coveragelist.size()-1);
	}


	/**
	 * Computes the number of extra lineages given a cluster (in a species tree) and a rooted gene tree when multiple alleles are sampled per species
	 *
	 * @param tr		a rooted gene tree
	 * @param cluster	the cluster being computed
	 * @param taxonMap	Maps gene tree taxa to species tree taxa
	 *
	 * @return the number of extra lineage
	 */
	public static int getClusterCoalNum_rooted(Tree tr, STITreeCluster cluster, Map<String, String> taxonMap) {
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		List<String> taxa = new LinkedList<String>();	// List of species taxa.

		for (String t : cluster.getTaxa()) {
			taxa.add(t);
		}

		int count = 0;
		for (TNode node : tr.postTraverse()) {
			if (node.isLeaf()) {
				String stTaxon = taxonMap.get(node.getName());	// Get the corresponding species name.
				int index = taxa.indexOf(stTaxon);
				BitSet bs = new BitSet(taxa.size());
				bs.set(index);
				if (cluster.containsCluster(bs)) {
					count++;
				}

				map.put(node, bs);
			}
			else {
				BitSet bs = new BitSet(taxa.size());
				int intersect = 0;
				int childCount = node.getChildCount();
				for (TNode child : node.getChildren()) {
					BitSet v = map.get(child);
					bs.or(v);
					if(childCount>2){
						if(cluster.containsCluster(v)){
							intersect ++;
						}
					}
				}

				if (cluster.containsCluster(bs)) {
					count -= node.getChildCount();
					count++;
				}
				else if(intersect>1){
					count -= intersect;
					count ++;
				}

				map.put(node, bs);
			}
		}
		return Math.max(count-1, 0);
	}


	/**
	 * Computes the number of extra lineages given a cluster (in a species tree) and an unrooted gene tree when multiple alleles are sampled per species
	 *
	 * @param tr		an unrooted gene tree
	 * @param cluster	the cluster being computed
	 * @param taxonMap	Maps gene tree taxa to species tree taxa
	 *
	 * @return the number of extra lineage
	 */
	public static int getClusterCoalNum_unrooted(Tree tr, STITreeCluster cluster, Map<String,String> taxonMap) {
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		List<String> gtTaxalist = new ArrayList<String>();
		String[] gtTaxa = tr.getLeaves();
		int ngtTaxa = gtTaxa.length;
		for(String leaf: gtTaxa){
			gtTaxalist.add(leaf);
		}
		STITreeCluster concluster = new STITreeCluster(gtTaxa);
		for(TNode n: tr.getNodes()){
			if(n.isLeaf()){
				if(cluster.containsLeaf(taxonMap.get(n.getName()))){
					concluster.addLeaf(n.getName());
				}
			}
		}

		List<BitSet> coveragelist = new ArrayList<BitSet>();
		for (int i = concluster.getCluster().nextSetBit(0); i >= 0; i = concluster.getCluster().nextSetBit(i+1)) {
			BitSet bs = new BitSet(ngtTaxa);
			bs.set(i);
			coveragelist.add(bs);
		}
		for (TNode node :tr.postTraverse()) {
			if(coveragelist.size()<=1){
				break;
			}
			BitSet bs = new BitSet(ngtTaxa);
			int intersect = 0;
			BitSet virtualbs = new BitSet(ngtTaxa);
			if (node.isLeaf()) {
				int index = gtTaxalist.indexOf(node.getName());
				bs.set(index);
				map.put(node, bs);
			}
			else {
				for (TNode child : node.getChildren()) {
					BitSet v = map.get(child);
					bs.or(v);
					if(concluster.containsCluster(v)){
						intersect ++;
						virtualbs.or(v);
					}
				}

				if (concluster.containsCluster(bs) || intersect>1) {
					if(concluster.containsCluster(bs)){
						virtualbs = bs;
					}
					for(int i=0; i<coveragelist.size();i++){
						BitSet exbs = coveragelist.get(i);
						BitSet temp = (BitSet) virtualbs.clone();
						temp.and(exbs);
						if(temp.equals(exbs)){
							coveragelist.remove(i);
							i--;
						}
					}
					coveragelist.add(virtualbs);
				}

				map.put(node, bs);
			}

			if(!node.isRoot()){
				BitSet complementbs = (BitSet)bs.clone();
				complementbs.flip(0, gtTaxa.length);
				if(intersect>0){
					complementbs.or(virtualbs);
				}
				if(concluster.containsCluster(complementbs)){
					for(int i=0; i<coveragelist.size();i++){
						BitSet exbs = coveragelist.get(i);
						BitSet temp = (BitSet) complementbs.clone();
						temp.and(exbs);
						if(temp.equals(exbs)){
							coveragelist.remove(i);
							i--;
						}
					}
					coveragelist.add(complementbs);
					break;
				}
			}
		}
		return Math.max(0,coveragelist.size()-1);
	}



}
