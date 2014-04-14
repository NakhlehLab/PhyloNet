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
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/26/11
 * Time: 4:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneTreeRefinement {
	/**
	 * Refine and root a set of gene trees optimally with respect to a given species tree
	 *
	 * @param gts	a list of given gene trees
	 * @param st	a given species tree
	 * @param taxonMap	Maps gene tree taxa to species tree taxa.
	 * @param rooted	whether treating the gene trees and species tree rooted or not
	 *
	 * @return the number of extra lineage
	 */
	public static void processGeneTrees(List<Tree> gts,Tree st, Map<String, List<String>> species2alleles, boolean rooted, double bootstrap){
		if(bootstrap<1){
			for(Tree tr: gts){
				if(Trees.handleBootStrapInTree(tr, bootstrap)==-1){
					throw new IllegalArgumentException("Input gene trees have nodes that don't have bootstrap value");
				}
			}
		}

		HashSet<STITreeCluster> speciesClusters = null;
		Tree st_copy = new STITree(st);

		if(species2alleles != null){
			List<String> taxa_list = new ArrayList<String>();
			for(Map.Entry<String, List<String>> entry: species2alleles.entrySet()){
				taxa_list.addAll(entry.getValue());
			}
			String[] taxa = new String[taxa_list.size()];
			taxa_list.toArray(taxa);
			taxa_list.clear();

			speciesClusters = new HashSet<STITreeCluster>();

			List<TNode> leafNodes = new ArrayList<TNode>();
			for(TNode node: st_copy.getNodes()){
				if(node.isLeaf()){
					leafNodes.add(node);
				}

			}
			for(TNode leaf: leafNodes){
				((STINode)leaf).setName("species_" + leaf.getName());
			}
			for(TNode leaf: leafNodes){
				String originalName = leaf.getName().substring(8);
				List<String> alleles = species2alleles.get(originalName);
				if(alleles == null){
					throw new IllegalArgumentException("The taxon " + leaf + " in the species tree is not in the mapping file");
				}
				if(alleles.size() == 1){
					((STINode)leaf).setName(alleles.get(0));
				}
				else{
					STITreeCluster c = new STITreeCluster(taxa);
					for(String allele: alleles){
						((STINode)leaf).createChild(allele);
						c.addLeaf(allele);
					}
					speciesClusters.add(c);
				}
			}
		}

		/*
		// for check correctness
		SymmetricDifference sd = new SymmetricDifference();
		int[] fn = new int[gts.size()];
		int[] fp = new int[gts.size()];
		int index = 0;
		for(Tree gt: gts){
			sd.computeDifference(gt, st_copy);
			fn[index] = sd.getFalseNegativeCount();
			fp[index] = computeFPUnRootedDistance(gt, st_copy, speciesClusters);
			//fp[index] = DeepCoalescencesCounter.countExtraCoal(gts, st_copy, rooted, bootstrap);
			index++;
		}
		*/

		if(!rooted){
			rootGeneTrees(gts, st_copy);
		}
		refineRootedGeneTrees(gts, st_copy, speciesClusters);

		/*
		// for check correctness
		index = 0;
		for(Tree gt: gts){
			//System.out.println(gt);
			int nfn = computeFNRootedDistance(gt, st_copy);
			int nfp = computeFPRootedDistance(gt, st_copy, speciesClusters);
			//int nfp = DeepCoalescencesCounter.countExtraCoal(gts, st_copy, true, bootstrap);
			if(fp[index]==0){
				if(nfn!=0 || nfp!=0){
					System.out.println(path);
					System.out.println("Error gene trees: " + gt);
				}
			}
			else{
				if(fp[index]!=nfp || nfn > fn[index]){
					System.out.println(path);
					System.out.println("Error gene trees: #" + (index+1) + " " + gt);
					System.out.println("fp:" + fp[index] + "/" + nfp);
					System.out.println("fn:" + fn[index] + "/" + nfn);
				}
			}
			index++;
		}
		*/
	}



	/**
	 * Refine a set of rooted gene trees with respect to a species tree
	 *
	 * @param gts		a set of gene trees
	 * @param st		a given species tree
	 *
	 */
	public static void refineRootedGeneTrees(List<Tree> gts, Tree st, HashSet<STITreeCluster> speciesClusters) {
		String[] st_taxa;
		if(speciesClusters == null){
			st_taxa = st.getLeaves();
		}
		else{
			st_taxa = speciesClusters.iterator().next().getTaxa();
		}
		List<STITreeCluster> original_clusters = st.getClusters(st_taxa, false);

		for(Tree gt: gts){
			Tree st_pruned = new STITree(st);
			List<STITreeCluster> st_clusters;
			String[] taxa = gt.getLeaves();
			if(taxa.length < st_taxa.length){
				((MutableTree)st_pruned).constrainByLeaves(Arrays.asList(taxa));
				st_clusters = st_pruned.getClusters(taxa, false);
			}
			else{
				taxa = st_taxa;
				st_clusters = original_clusters;
			}
			if(!Trees.leafSetsAgree(gt, st_pruned)){
				throw new RuntimeException("Gene tree " + gt.toNewick() + " has some taxa that species tree doesn't have");
			}
			List<List<TNode>> refineNodes = new ArrayList<List<TNode>>();

			Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();

			for (TNode node : gt.postTraverse()) {
				if (node.isLeaf()) {
					int index = -1;
					for(int i=0; i<taxa.length; i++){
						if(taxa[i].equals(node.getName())){
							index = i;
							break;
						}
					}
					BitSet bs = new BitSet(taxa.length);
					bs.set(index);
					map.put(node, bs);
				}
				else {
					BitSet bs = new BitSet();
					int childCount = 0;
					for (TNode child : node.getChildren()) {
						BitSet v = map.get(child);
						bs.or(v);
						childCount ++;
					}
					if(childCount > 2){
						STITreeCluster bsCluster = new STITreeCluster(st_taxa);
						bsCluster.setCluster(bs);

						for(STITreeCluster cluster: st_clusters){
							boolean isSpeciesCluster = false;
							if(speciesClusters!=null){
								isSpeciesCluster = speciesClusters.contains(cluster);
							}

							if((cluster.containsCluster(bs)&&!isSpeciesCluster) || cluster.isDisjoint(bs)){
								continue;
							}
							if(!bsCluster.containsCluster(cluster) && !isSpeciesCluster){
								continue;
							}

							List<TNode> refineCandidate = new ArrayList<TNode>();

							BitSet union = new BitSet(taxa.length);
							for (TNode child : node.getChildren()) {
								BitSet v = map.get(child);
								if(cluster.containsCluster(v)){
									union.or(v);
									refineCandidate.add(child);
								}
							}

							if(refineCandidate.size() > 1 && ( union.equals(cluster.getCluster()) || isSpeciesCluster)){
								refineNodes.add(refineCandidate);
							}
						}
					}
					map.put(node, bs);
				}
			}


			for(int i=1;i<refineNodes.size();i++){
				List<TNode> e1 = refineNodes.get(i);
				for(int j=0;j<i;j++){
					List<TNode> e2 = refineNodes.get(j);
					if(e1.size() > e2.size()){
						refineNodes.remove(i);
						refineNodes.add(j,e1);
						break;
					}
				}
			}

			for(List<TNode> toRefine: refineNodes){
				TNode parent = toRefine.get(0).getParent();
				STINode new_parent = ((STINode)parent).createChild();
				for(TNode node: toRefine){
					new_parent.adoptChild((TMutableNode)node);
				}
			}

			Trees.removeBinaryNodes((MutableTree)gt);
		}
	}




	/**
	 * Root a set of unrooted gene trees with respect to a species tree
	 *
	 * @param gts		a set of gene trees
	 * @param st		a given species tree
	 *
	 */
	public static void rootGeneTrees(List<Tree> gts, Tree st) {
		String[] st_taxa = st.getLeaves();
		List<STITreeCluster> original_clusters = st.getClusters(st_taxa, false);
		for(int i=1;i<original_clusters.size();i++){
			STITreeCluster e1 = original_clusters.get(i);
			int size1 = e1.getClusterSize();
			for(int j=0;j<i;j++){
				STITreeCluster e2 = original_clusters.get(j);
				if( size1 > e2.getClusterSize()){
					original_clusters.remove(i);
					original_clusters.add(j,e1);
					break;
				}
			}
		}

		for(Tree gt: gts){
			Tree st_pruned = new STITree(st);
			List<STITreeCluster> st_clusters;
			String[] taxa = gt.getLeaves();
			if(taxa.length < st_taxa.length){
				((MutableTree)st_pruned).constrainByLeaves(Arrays.asList(taxa));
				st_clusters = st_pruned.getClusters(taxa, false);
				for(int i=1;i<st_clusters.size();i++){
					STITreeCluster e1 = st_clusters.get(i);
					int size1 = e1.getClusterSize();
					for(int j=0;j<i;j++){
						STITreeCluster e2 = st_clusters.get(j);
						if( size1 > e2.getClusterSize()){
							st_clusters.remove(i);
							st_clusters.add(j,e1);
							break;
						}
					}
				}
			}
			else{
				taxa = st_taxa;
				st_clusters = original_clusters;
			}
			if(!Trees.leafSetsAgree(gt, st_pruned)){
				throw new RuntimeException("Gene tree " + gt.toNewick() + " has some taxa that species tree doesn't have");
			}

			if(st_clusters.get(0).getClusterSize() == taxa.length-1){
				String leaf = taxa[st_clusters.get(0).getCluster().nextClearBit(0)];
				gt.rerootTreeAtEdge(leaf);
				continue;
			}

			List<STITreeCluster> gt_clusters = gt.getBipartitionClusters(taxa, false);
			for(int i=1;i<gt_clusters.size();i++){
				STITreeCluster e1 = gt_clusters.get(i);
				int size1 = e1.getClusterSize();
				for(int j=0;j<i;j++){
					STITreeCluster e2 = gt_clusters.get(j);
					if( size1 > e2.getClusterSize()){
						gt_clusters.remove(i);
						gt_clusters.add(j,e1);
						break;
					}
				}
			}
			STITreeCluster max_cluster = null;
			for(STITreeCluster gtcl: gt_clusters){
				if(gtcl.getClusterSize() == 1){
					break;
				}
				int gsize = gtcl.getClusterSize();
				for(STITreeCluster stcl: st_clusters){
					if(stcl.getClusterSize() < gsize){
						break;
					}
					if(stcl.containsCluster(gtcl)){
						max_cluster = gtcl;
					}
				}
				if(max_cluster != null){
					break;
				}
			}

			if(max_cluster == null){
				if(gt.getRoot().getChildCount() == 2){
					for(TNode root_child: gt.getRoot().getChildren()){
						if(!root_child.isLeaf()){
							gt.rerootTreeAtNode(root_child);
							break;
						}
					}
				}
				continue;
			}

			STITreeCluster cmax_cluster = max_cluster.complementaryCluster();
			TNode max_node = null;
			Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
			for (TNode node : gt.postTraverse()) {
				if (node.isLeaf()) {
					int index = -1;
					for(int i=0; i<taxa.length; i++){
						if(taxa[i].equals(node.getName())){
							index = i;
							break;
						}
					}
					BitSet bs = new BitSet(taxa.length);
					bs.set(index);
					map.put(node, bs);
				}
				else {
					BitSet bs = new BitSet(taxa.length);
					for (TNode child : node.getChildren()) {
						BitSet v = map.get(child);
						bs.or(v);
					}
					if(bs.equals(max_cluster.getCluster())){
						max_node = node.getParent();
						break;
					}
					else{
						if(bs.equals(cmax_cluster.getCluster())){
							max_node = node;
							break;
						}
					}
					map.put(node, bs);
				}
			}
			gt.rerootTreeAtNode(max_node);
		}
	}







}
