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
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by Yun Yu
 * Date: 11/9/11
 * Time: 4:57 PM
 *
 * This class is to infer species trees from a collection of rooted gene trees with branch lengths under MDC
 *
 * See "From Gene Trees to Species Trees: Algorithms for Parsimonious Reconciliation", Yun Yu's master thesis 2012
 */
public class MDCWithTime {

	/**
	 * Infers the species tree from a collection of unltrametric rooted gene trees with single allele sampled per species by dynamic programming.
	 *
	 * @param 	trees	the collection of ultrametric gene trees (with branch lengths) from which the species tree is to be inferred
     *
	 * @return	the tree inferred from the supplied gene trees with the total number of extra lineages
	 */
	public Solution inferSpeciesTree(List<MutableTuple<Tree,Double>> trees){
		if (trees == null || trees.size() == 0) {
			System.err.println("Empty list of trees. The function returns a null tree.");
			return null;
		}

		List<String> taxalist = new ArrayList<String>();
		for(MutableTuple<Tree,Double> tr: trees){
			for (TNode node : tr.Item1.postTraverse()) {
				if (node.isLeaf() && taxalist.indexOf(node.getName())==-1) {
					taxalist.add(node.getName());
				}
			}
            setNodeData(tr.Item1);
		}

		String[] taxa = new String[taxalist.size()];
		int index = 0;
		for(String taxon: taxalist){
			taxa[index++] = taxon;
		}

		Map<Integer, List<Branch>> branches = new HashMap<Integer, List<Branch>>();
		double maxEL = computeBranchInGraph(trees,taxa,branches);
		return findTreesByDP(branches,taxa,maxEL);
	}


	/**
	 * Infers the species tree from a collection of unltrametric rooted gene trees with multiple alleles sampled per species by dynamic programming.
	 *
	 * @param 	trees	the collection of ultrametric gene trees (with branch lengths) from which the species tree is to be inferred
	 * @param	taxonMap	maps gene tree taxa to species tree taxa.
     *
	 * @return	the tree inferred from the supplied gene trees with the total number of extra lineages
	 *
	 */
	public Solution inferSpeciesTree(List<MutableTuple<Tree,Double>> trees,Map<String,String> taxonMap){
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


		for(MutableTuple<Tree,Double> tr:trees){
			setNodeData(tr.Item1);
		}
		Map<Integer, List<Branch>> branches = new HashMap<Integer, List<Branch>>();
		double maxEL = computeBranchInGraph(trees,stTaxa,gtTaxa,taxonMap,branches);

		return findTreesByDP(branches,stTaxa,maxEL);
	}


	/**
	 * Finds the species tree by dynamic programming
	 *
	 * @param 	branches	candidate branches with extra lineages that will be used in the species tree
	 * @param	stTaxa	  taxa expression
	 * @param   maxEL	  the maximal number of extra lineage
     *
	 * @return	the tree inferred from the supplied gene trees with the total number of extra lineages
	 *
	 */
	private Solution findTreesByDP(Map<Integer, List<Branch>> branches, String[] stTaxa, double maxEL){
		Branch all = branches.get(stTaxa.length).get(0);
		computeMinCost(all,branches,maxEL);

		// Build the minimum tree.
		List<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		List<Double> coals = new ArrayList<Double>();
		Stack<Branch> minBranches = new Stack<Branch>();

		if (all.getMinRightBranch() != null) {
			minBranches.push(all.getMinRightBranch());
		}
		if (all.getMinLeftBranch() != null) {
			minBranches.push(all.getMinLeftBranch());
		}
		if (all._subBrs != null) {
			for(Branch b: all._subBrs){
				minBranches.push(b);
			}
		}

		while (!minBranches.isEmpty()) {
			Branch pe = minBranches.pop();
			minClusters.add(pe.getChildCluster());
			coals.add(pe.getExtraLineage());
			if (pe.getMinRightBranch() != null) {
				minBranches.push(pe.getMinRightBranch());
			}
			if (pe.getMinLeftBranch() != null) {
				minBranches.push(pe.getMinLeftBranch());
			}
			if (pe._subBrs != null) {
				for(Branch b: pe._subBrs){
					minBranches.push(b);
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
			sol._st =  Trees.buildTreeFromClusters(minClusters);
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

		sol._totalCoals=all._min_cost;
		return sol;
	}


	/**
	 * Compute the minimum cost for trees whose leaf set is branches's child cluster.
	 *
	 * @param	branches	candidate branches with extra lineages that will be used in the species tree
	 * @param	br	the branch being computed
	 * @param	maxEL  the maximal number of extra lineages
     *
	 * @return  maximal number of extra lineages
	 */
	private double computeMinCost(Branch br,Map<Integer, List<Branch>> branches, double maxEL){
		if (br.getMaxScore() != -1) {
			return br.getMaxScore();
		}

		STITreeCluster<Double> cl = br.getChildCluster();
		int vsize = cl.getClusterSize();

		if (vsize <= 1) {
			br.setMinCost(br.getExtraLineage());
			br.setMaxScore(maxEL - br.getExtraLineage());
			br.setMinLeftBranch(null);
			br.setMinRightBranch(null);

			return br.getMaxScore();
		}

		for (int i = 1; i <= (vsize / 2); i++) {
			List<Branch> leftList = branches.get(i);
			if (leftList != null) {
				for (Branch lb : leftList) {
					if (cl.containsCluster(lb.getChildCluster())
							&& lb.getParentBitSet().equals(cl.getCluster())
							&& lb.getParentCluster().getData()<=cl.getData()) {
						double lscore = computeMinCost(lb,branches,maxEL);
						double lcost = lb.getMinCost();

						List<Branch> rightList = branches.get(vsize - i);
						if (rightList != null) {
							for (Branch rb : rightList) {
								if (lb.getChildCluster().isDisjoint(rb.getChildCluster()) && cl.containsCluster(rb.getChildCluster())
										&& rb.getParentBitSet().equals(cl.getCluster()) && rb.getParentCluster().getData()<=cl.getData()) {
									double rscore = computeMinCost(rb,branches,maxEL);
									double rcost = rb.getMinCost();
									if ((br.getMaxScore() == -1) || (lscore + rscore + maxEL - br.getExtraLineage() > br.getMaxScore())) {
										br.setMaxScore(lscore + rscore + maxEL - br.getExtraLineage());
										br.setMinCost(lcost + rcost + br.getExtraLineage());
										br.setMinLeftBranch(lb);
										br.setMinRightBranch(rb);
									}
									break;	// Already found the only pair of clusters whose union is v's cluster.
								}
							}
						}
					}
				}
			}
		}

		if (br._max_score == -1) {	// Cannot resolve v as a binary node.

			ArrayList<Branch> subBranches = new ArrayList<Branch>();
			//get all the subbranches
			for(int i=1; i<cl.getClusterSize(); i++){
				List<Branch> brs = branches.get(i);
				if(brs == null)continue;
				for(Branch tb: brs){
					if(cl.containsCluster(tb.getParentCluster()) && tb.getParentCluster().getData()<=cl.getData()){
						if(tb.getMaxScore()==-1){
							computeMinCost(tb, branches, maxEL);
						}
						subBranches.add(tb);
					}
				}
			}


			for(int i=1;i<subBranches.size();i++){
				Branch b1 = subBranches.get(i);
				for(int j=0;j<i;j++){
					Branch b2 = subBranches.get(j);
					if(b1.getChildCluster().getClusterSize() == b2.getChildCluster().getClusterSize()){
						if(b1.getParentCluster().getClusterSize() == b2.getParentCluster().getClusterSize()){
							if(b1._max_score > b2._max_score){
								subBranches.remove(i);
								subBranches.add(j,b1);
								break;
							}
							else if(b1._max_score == b2._max_score){
								if(b1.equals(b2) && b1.getParentCluster().getData()< b2.getParentCluster().getData()){
									subBranches.remove(i);
									subBranches.add(j,b1);
									break;
								}
							}
						}
						else if(b1.getParentCluster().getClusterSize() < b2.getParentCluster().getClusterSize()){
							subBranches.remove(i);
							subBranches.add(j,b1);
							break;
						}
					}
					else if(b1.getChildCluster().getClusterSize() > b2.getChildCluster().getClusterSize()){
						subBranches.remove(i);
						subBranches.add(j,b1);
						break;
					}
				}
			}

			ArrayList<Branch> maxSubBranches = new ArrayList<Branch>();
			for(Branch tb: subBranches){
				boolean compatible = true;
				boolean contained = false;
				boolean exist = false;
				ArrayList<Branch> rmBrs = new ArrayList<Branch>();
				for(Branch ex_b: maxSubBranches){
					if(tb.equals(ex_b)){
						exist = true;
						break;
					}
					if(ex_b.containsBranch(tb)){
						contained = true;
						break;
					}
					if(!tb.isCompatible(ex_b)){
						compatible = false;
						break;
					}
					if(tb.containsBranch(ex_b)){
						rmBrs.add(ex_b);

					}
				}
				if(compatible && (!contained) && (!exist)){
					for(Branch ex_b: rmBrs){
						maxSubBranches.remove(ex_b);
					}
					maxSubBranches.add(tb);
				}
			}

			br._min_cost = br.getExtraLineage();
			br._max_score = maxEL - br.getExtraLineage();
			for(Branch bv: maxSubBranches){
				br._min_cost += bv._min_cost;
				br._max_score += bv._max_score;
			}
			br._subBrs = maxSubBranches;
		}

		return br.getMaxScore();
	}


	/**
	 * Constructs the branches with time and their number of extra lineages when single allele is sampled per species
	 *
	 * @param	trees	a list of given gene trees (ultrametric)
	 * @param	taxa	taxa in species/gene trees
	 * @param	branches	resulting branches with extra lineages
	 *
	 * @return	maximal number of extra lineages
	 */
	private double computeBranchInGraph(List<MutableTuple<Tree,Double>> trees,String[] taxa, Map<Integer, List<Branch>> branches){
		List<STITreeCluster<Double>> clusters1 = computeTreeClustersWithTime(trees,taxa);
		HashMap<STITreeCluster,Double> taxonPairTime = computeTaxonPairTime(clusters1,taxa);

		LinkedList<STITreeCluster<Double>> clusters2 = new LinkedList<STITreeCluster<Double>>();
		for(STITreeCluster<Double> cl:clusters1){
			if(cl.getClusterSize()>1)
				clusters2.add(cl);
		}

		List<Branch> l = new LinkedList<Branch>();
		Branch root_branch = new Branch(clusters1.get(0),null);
		root_branch.setExtraLineage(0);
		root_branch.setMinCost(-1);
		l.add(root_branch);
		branches.put(taxa.length, l);

		double maxEL = 0;
		for(STITreeCluster<Double> c1:clusters1){
			int c1_size = c1.getClusterSize();
			l = branches.get(c1_size);
			if(l == null){
				l = new LinkedList<Branch>();
				branches.put(c1_size, l);
			}

			for(STITreeCluster<Double> c2:clusters2){
				int c2_size = c2.getClusterSize();
				if(c2_size == c1_size){
					break;
				}
				if(c2.containsCluster(c1)){
				//if(c1.canMakeBranch(c2, clusters2)){
					List<Branch> ex_l = branches.get(c2_size);
					List<Double> time_l = new LinkedList<Double>();
					for(Branch b:ex_l){
						STITreeCluster<Double> ex_c1 = b.getChildCluster();
						if(ex_c1.getCluster().equals(c2.getCluster())
								&& (!time_l.contains(ex_c1.getData()))){
							time_l.add(ex_c1.getData());
						}
					}
					if(time_l.size()==0){
						Branch newBranch = new Branch(c1.duplicate(),c2.duplicate());
						reviseTime(newBranch,taxonPairTime);
						if(c1_size>1){
							double el = getBranchCoalNum(newBranch,trees);
							newBranch.setExtraLineage(el);
							if(el > maxEL){
								maxEL = el;
							}
						}
						else{
							newBranch.setExtraLineage(0);
						}
						newBranch.setMinCost(-1);
						//l.add(newBranch);
						int pos = l.indexOf(newBranch);
						if(pos == -1){
							l.add(newBranch);
						}else{
							double time = newBranch.getParentCluster().getData();
							double el = newBranch.getExtraLineage();
							while(pos != l.size()){
								Branch ex_b = l.get(pos);
								if(!ex_b.equals(newBranch)){
									l.add(pos, newBranch);
									break;
								}
								else if((ex_b.getExtraLineage() > el) || (ex_b.getExtraLineage() == el && ex_b.getParentCluster().getData() > time)){
									l.add(pos, newBranch);
									break;
								}
								else{
									pos++;
								}
							}
							if(pos == l.size()){
								l.add(newBranch);
							}
						}

					}
					else{
						List<Double> tl = new LinkedList<Double>();
						for(double t:time_l){
							Branch newBranch = new Branch(c1.duplicate(),c2.duplicate());
							newBranch.getParentCluster().setData(t);
							reviseTime(newBranch,taxonPairTime);
							if(tl.contains(newBranch.getParentCluster().getData())){
								continue;
							}
							else{
								if(c1_size>1){
									double el = getBranchCoalNum(newBranch,trees);
									newBranch.setExtraLineage(el);
									if(el > maxEL){
										maxEL = el;
									}
								}
								else{
									newBranch.setExtraLineage(0);
								}
								newBranch.setMinCost(-1);
								int pos = l.indexOf(newBranch);
								if(pos == -1){
									l.add(newBranch);
								}else{
									double time = newBranch.getParentCluster().getData();
									double el = newBranch.getExtraLineage();
									while(pos != l.size()){
										Branch ex_b = l.get(pos);
										if(!ex_b.equals(newBranch)){
											l.add(pos, newBranch);
											break;
										}
										else if((ex_b.getExtraLineage() > el) || (ex_b.getExtraLineage() == el && ex_b.getParentCluster().getData() > time)){
											l.add(pos, newBranch);
											break;
										}
										else{
											pos++;
										}
									}
									if(pos == l.size()){
										l.add(newBranch);
									}
								}
								tl.add(newBranch.getParentCluster().getData());
							}
						}
					}
				}
			}
		}
		maxEL = maxEL+1;
		return maxEL;
	}


	/**
	 * Constructs the branches with time and their number of extra lineages when multiple alleles are sampled per species
	 *
	 * @param	trees	a list of given gene trees (ultrametric)
	 * @param	stTaxa	taxa in species tree
	 * @param	gtTaxa	taxa in gene trees
	 * @param	taxonMap	maps gene tree taxa to species tree taxa
	 * @param	branches	resulting branches with extra lineages
	 *
	 * @return	maximal number of extra lineages
	 */
	private double computeBranchInGraph(List<MutableTuple<Tree,Double>> trees,String stTaxa[],String gtTaxa[],Map<String,String> taxonMap, Map<Integer, List<Branch>> branches){
		List<STITreeCluster<Double>> clusters1 = computeTreeClustersWithTime(trees,stTaxa,gtTaxa,taxonMap);
		HashMap<STITreeCluster,Double> taxonPairTime = computeTaxonPairTime(clusters1,stTaxa);

		LinkedList<STITreeCluster<Double>> clusters2 = new LinkedList<STITreeCluster<Double>>();
		for(STITreeCluster<Double> cl:clusters1){
			if(cl.getClusterSize()>1)
				clusters2.add(cl);
		}

		List<Branch> l = new LinkedList<Branch>();
		Branch root_branch = new Branch(clusters1.get(0),null);
		root_branch.setExtraLineage(0);
		root_branch.setMinCost(-1);
		l.add(root_branch);
		branches.put(stTaxa.length, l);
		double maxEL = 0;
		LinkedList<STITreeCluster<Double>> uselessClusters = new LinkedList<STITreeCluster<Double>>();
		for(STITreeCluster<Double> c1:clusters1){
			boolean added = false;
			int c1_size = c1.getClusterSize();
			l = branches.get(c1_size);
			if(l == null){
				l = new LinkedList<Branch>();
				branches.put(c1_size, l);
			}

			for(STITreeCluster<Double> c2:clusters2){
				int c2_size = c2.getClusterSize();
				if(c2_size == c1_size){
					break;
				}
				if(uselessClusters.contains(c2))continue;
				if(c2.containsCluster(c1)){
				//if(c1.canMakeBranch(c2, clusters2)){
					List<Branch> ex_l = branches.get(c2_size);
					List<Double> time_l = new LinkedList<Double>();
					for(Branch b:ex_l){
						STITreeCluster<Double> ex_c1 = b.getChildCluster();
						if(ex_c1.getCluster().equals(c2.getCluster())
								&& (!time_l.contains(ex_c1.getData()))){
							time_l.add(ex_c1.getData());
						}
					}
					if(time_l.size()==0){
						Branch newBranch = new Branch(c1.duplicate(),c2.duplicate());
						reviseTime(newBranch,taxonPairTime);
						double el = getBranchCoalNum(newBranch,trees,taxonMap);
						newBranch.setExtraLineage(el);
						if(el > maxEL){
							maxEL = el;
						}
						newBranch.setMinCost(-1);
						int pos = l.indexOf(newBranch);
						if(pos == -1){
							l.add(newBranch);
						}else{
							double time = newBranch.getParentCluster().getData();
							while(pos != l.size()){
								Branch ex_b = l.get(pos);
								if(!ex_b.equals(newBranch)){
									l.add(pos, newBranch);
									break;
								}
								else if((ex_b.getExtraLineage() > el) || (ex_b.getExtraLineage() == el && ex_b.getParentCluster().getData() > time)){
									l.add(pos, newBranch);
									break;
								}
								else{
									pos++;
								}
							}
							if(pos == l.size()){
								l.add(newBranch);
							}
						}
						added=true;
					}
					else{
						List<Double> tl = new LinkedList<Double>();
						for(double t:time_l){
							Branch newBranch = new Branch(c1.duplicate(),c2.duplicate());
							newBranch.getParentCluster().setData(t);
							reviseTime(newBranch,taxonPairTime);
							if(tl.contains(newBranch.getParentCluster().getData())){
								continue;
							}
							else{
								double el = getBranchCoalNum(newBranch,trees,taxonMap);
								newBranch.setExtraLineage(el);
								if(el > maxEL){
									maxEL = el;
								}
								newBranch.setMinCost(-1);
								//l.add(newBranch);
								int pos = l.indexOf(newBranch);
								if(pos == -1){
									l.add(newBranch);
								}else{
									double time = newBranch.getParentCluster().getData();
									while(pos != l.size()){
										Branch ex_b = l.get(pos);
										if(!ex_b.equals(newBranch)){
											l.add(pos, newBranch);
											break;
										}
										else if((ex_b.getExtraLineage() > el) || (ex_b.getExtraLineage() == el && ex_b.getParentCluster().getData() > time)){
											l.add(pos, newBranch);
											break;
										}
										else{
											pos++;
										}
									}
									if(pos == l.size()){
										l.add(newBranch);
									}
								}
								added=true;
								tl.add(newBranch.getParentCluster().getData());
							}
						}
					}
				}
			}
			if(added==false){
				if(c1.getClusterSize()!=stTaxa.length)
					uselessClusters.add(c1);
			}
		}
		maxEL = maxEL+1;
		return maxEL;
	}


	/**
	 * Revises the time of a branch (the height of the node that this branch is incident with) according to constraints
	 *
	 * @param	br	the branch being revised
	 * @param	taxonPairTime	minimal pairwise coalescence time
	 *
	 */
	private void reviseTime(Branch br, HashMap<STITreeCluster,Double> taxonPairTime){
		double epsilon =0.0;
		STITreeCluster<Double> child = br.getChildCluster();
		STITreeCluster<Double> parent = br.getParentCluster();
		String[] leaves1 = child.getClusterLeaves();
		String[] leaves2 = new String[parent.getClusterSize()-child.getClusterSize()];
		BitSet bs2 = (BitSet)(parent.getCluster().clone());
		bs2.xor(child.getCluster());
		String[] taxa = child.getTaxa();
		int index = 0;
		for(int i=0;i<taxa.length;i++){
			if(bs2.get(i)){
				leaves2[index++] = taxa[i];
			}
		}


		double minTime = -1;
		for(String taxonA: leaves1){
			for(String taxonB: leaves2){
				STITreeCluster cl = new STITreeCluster(br.getTaxa());
				cl.addLeaf(taxonA);
				cl.addLeaf(taxonB);

				double time = taxonPairTime.get(cl);
				if(minTime==-1 || time<minTime)
					minTime = time;
			}
		}
		if(parent.getData()>minTime){
			parent.setData(minTime);
		}
		if(child.getData()>=parent.getData()){
			child.setData(parent.getData()-epsilon);
		}
	}

	/**
	 * Computes the number of extra lineages contributed by a branch in a species tree and a collection of gene trees when single allele is sampled per species
	 *
	 * @param br	the branch being computed
	 * @param trees	a list of given gene trees
	 *
	 * @return the number of extra lineage
	 */
	private double getBranchCoalNum(Branch br,List<MutableTuple<Tree,Double>> trees){
		double coalNum = 0;

		if(br.getChildCluster().getClusterSize()==1
				&& br.getParentCluster().getClusterSize()==2){
			coalNum = 0;
		}
		else{
			for (MutableTuple<Tree,Double> tr : trees) {
				coalNum += getBranchCoalNum(br,tr.Item1) * tr.Item2;
			}
		}
		return coalNum;
	}


	/**
     * Computes the number of extra lineages contributed by a branch in a species tree and a gene trees when single allele is sampled per species
	 *
	 * @param br	the branch being computed
	 * @param tr	a given gene tree
	 *
	 * @return the number of extra lineage
	 */
	private int getBranchCoalNum(Branch br,Tree tr){
		int coalNum = 0;
		STITreeCluster<Double> child = br.getChildCluster();
		STITreeCluster<Double> parent =  br.getParentCluster();

		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		List<String> taxa = new LinkedList<String>();

		for (String t : child.getTaxa()) {
			taxa.add(t);
		}

		for (TNode node : tr.postTraverse()) {
			if (node.isLeaf()) {
				int index = taxa.indexOf(node.getName());
				BitSet bs = new BitSet();

				bs.set(index);
				if (child.containsCluster(bs)) {
					coalNum++;
				}

				map.put(node, bs);
			}
			else {
				BitSet bs = new BitSet(taxa.size());
				for (TNode ch : node.getChildren()) {
					BitSet v = map.get(ch);
					bs.or(v);
				}
				double coalTime = ((STINode<Double>)node).getData();
				if (child.containsCluster(bs) && coalTime < parent.getData()) {
					coalNum -= node.getChildCount();
					coalNum++;
				}

				map.put(node, bs);
			}
		}

		return Math.max(0, coalNum-1);
	}


	/**
     * Computes the number of extra lineages contributed by a branch in a species tree and a collection of gene trees when multiple alleles are sampled per species
	 *
	 * @param br	the branch being computed
	 * @param trees		a list of given gene trees
	 * @param taxonMap	maps gene tree taxa to species tree taxa.
	 *
	 * @return the number of extra lineage
	 */
	private double getBranchCoalNum(Branch br,List<MutableTuple<Tree,Double>> trees,Map<String,String> taxonMap){
		double coalNum = 0;

		for (MutableTuple<Tree,Double> tr : trees) {
			coalNum += getBranchCoalNum(br,tr.Item1,taxonMap) * tr.Item2;
		}

		return coalNum;
	}


	/**
     * Computes the number of extra lineages contributed by a branch in a species tree and a collection of gene trees when multiple alleles are sampled per species
	 *
	 * @param br	the branch being computed
	 * @param tr	a given gene tree
	 * @param taxonMap	maps gene tree taxa to species tree taxa.
	 *
	 * @return the number of extra lineage
	 */
	private int getBranchCoalNum(Branch br,Tree tr,Map<String,String> taxonMap){
		int coalNum = 0;
		STITreeCluster<Double> child = br.getChildCluster();
		STITreeCluster<Double> parent =  br.getParentCluster();

		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		List<String> taxa = new LinkedList<String>();

		for (String t : child.getTaxa()) {
			taxa.add(t);
		}

		for (TNode node : tr.postTraverse()) {
			if (node.isLeaf()) {
				String stTaxon = taxonMap.get(node.getName());
				int index = taxa.indexOf(stTaxon);
				BitSet bs = new BitSet();

				bs.set(index);
				if (child.containsCluster(bs)) {
					coalNum++;
				}

				map.put(node, bs);
			}
			else {
				BitSet bs = new BitSet(taxa.size());
				for (TNode ch : node.getChildren()) {
					BitSet v = map.get(ch);
					bs.or(v);
				}
				double coalTime = ((STINode<Double>)node).getData();
				if (child.containsCluster(bs) && coalTime < parent.getData()) {
					coalNum -= node.getChildCount();
					coalNum++;
				}

				map.put(node, bs);
			}
		}

		return Math.max(0, coalNum-1);
	}


	/**
	 * Computes the minimum pairwise distance of taxon pairs
	 *
	 * @param clusters	a list of all clusters with time
	 * @param stTaxa	taxa in species tree
	 *
	 * @return  maps from taxon pairs to distance
	 */
	private HashMap<STITreeCluster,Double> computeTaxonPairTime(List<STITreeCluster<Double>> clusters, String[] stTaxa){
		HashMap<STITreeCluster,Double> taxonPairTime = new HashMap<STITreeCluster, Double>();
		ArrayList<STITreeCluster<Double>> taxonPairs = new ArrayList<STITreeCluster<Double>>();

		for(int i=0;i<stTaxa.length;i++){
			for(int j=i+1;j<stTaxa.length;j++){
				STITreeCluster<Double> cl = new STITreeCluster<Double>(stTaxa);
				BitSet bs = new BitSet(stTaxa.length);
				bs.set(i);
				bs.set(j);
				cl.setCluster(bs);
				cl.setData(-1.0);
				taxonPairs.add(cl);
			}
		}

		for(STITreeCluster<Double> cl1: clusters){
			if(cl1.getClusterSize()==1)break;

			for(STITreeCluster<Double> cl2 : taxonPairs){
				if(cl1.containsCluster(cl2)){
					if(cl1.getData()<cl2.getData() || cl2.getData()==-1){
						cl2.setData(cl1.getData());
					}
				}
			}
		}

		for(STITreeCluster<Double> clwt : taxonPairs){
			STITreeCluster cl = new STITreeCluster(stTaxa);
			cl.setCluster(clwt.getCluster());
			taxonPairTime.put(cl, clwt.getData());
		}

		return taxonPairTime;
	}



	/**
	 * Computes clusters induced by gene trees when single allele is sampled per species
	 *
	 * @param trees	the list of gene trees
	 * @param taxa	taxa in species tree/gene trees
	 *
	 * @return  clusters with the number of extra lineages
	 */
	private List<STITreeCluster<Double>> computeTreeClustersWithTime(List<MutableTuple<Tree,Double>> trees,String taxa[]) {
		LinkedList<STITreeCluster<Double>> clusters = new LinkedList<STITreeCluster<Double>>();

		STITreeCluster<Double> all = new STITreeCluster<Double>(taxa);
		for (String t : taxa) {
			all.addLeaf(t);
		}
		double minTime = -1;
		for(MutableTuple<Tree,Double> tr:trees){
			TNode root = tr.Item1.getRoot();
			if(((STINode<Double>)root).getData() < minTime || minTime == -1){
				minTime = ((STINode<Double>)root).getData();
			}
		}
		all.setData(minTime);

		clusters.add(all);


        for(MutableTuple<Tree,Double> tr:trees){
			for (STITreeCluster<Double> tc : tr.Item1.getClusters(taxa, true)) {
				if(!clusters.contains(tc)){
					clusters.add(tc.duplicate());
				}
				else{
					int index = clusters.indexOf(tc);
					if(clusters.get(index).getData() > tc.getData()){
						clusters.get(index).setData(tc.getData());
					}
				}

			}
		}

		//order the clusters by size decreasing
		for(int i=1;i<clusters.size();i++){
			STITreeCluster<Double> cl1 = clusters.get(i);
			for(int j=0;j<i;j++){
				STITreeCluster<Double> cl2 = clusters.get(j);
				if(cl2.getClusterSize() < cl1.getClusterSize()){
					clusters.remove(cl1);
					clusters.add(j,cl1);
					break;
				}
			}
		}

		//adjust the time
		for(int i=0;i<clusters.size();i++){
			STITreeCluster<Double> cl1 = clusters.get(i);
			for(int j=i+1;j<clusters.size();j++){
				STITreeCluster<Double> cl2 = clusters.get(j);
				if(cl1.containsCluster(cl2)){
					if(cl1.getData() < cl2.getData()){
						cl2.setData(cl1.getData());
					}
				}
			}
		}

		 //Add clusters of size 1.
		for (String t : taxa) {
			STITreeCluster<Double> c = new STITreeCluster<Double>(taxa);
			c.addLeaf(t);
			c.setData(0.0);
			clusters.add(c);
		}


		return clusters;
	}


	/**
     * Computes clusters induced by gene trees when multiple alleles are sampled per species
	 *
	 * @param trees	the list of gene trees
	 * @param stTaxa	taxa in species tree
	 * @param gtTaxa	taxa in gene trees
	 * @param taxonMap	maps gene tree taxa to species tree taxa.
	 *
	 * @return	clusters with the number of extra lineages
	 */
	private List<STITreeCluster<Double>> computeTreeClustersWithTime(List<MutableTuple<Tree,Double>> trees,String stTaxa[],String gtTaxa[],Map<String,String> taxonMap) {
		LinkedList<STITreeCluster<Double>> clusters = new LinkedList<STITreeCluster<Double>>();

		//Add the cluster containing all taxa to the end of the list.
		STITreeCluster<Double> all = new STITreeCluster<Double>(stTaxa);
		for (String t : stTaxa) {
			all.addLeaf(t);
		}
		double minTime = -1;
        for(MutableTuple<Tree,Double> tr:trees){
			TNode root = tr.Item1.getRoot();
			if(((STINode<Double>)root).getData() < minTime || minTime == -1){
				minTime = ((STINode<Double>)root).getData();
			}
		}
		all.setData(minTime);

		clusters.add(all);

		// Compute the list of branches in the graph, which correspond to the induced clusters.
        for(MutableTuple<Tree,Double> tr:trees){
			for (STITreeCluster<Double> tc : tr.Item1.getClusters(gtTaxa, true)) {
				STITreeCluster<Double> stCluster = new STITreeCluster<Double>(stTaxa);
				for (String s : tc.getClusterLeaves()) {
					stCluster.addLeaf(taxonMap.get(s));
				}
				stCluster.setData(tc.getData());

				if(!clusters.contains(stCluster)){
					if(stCluster.getClusterSize()>1)
						clusters.add(stCluster);
				}
				else{
					int index = clusters.indexOf(stCluster);
					if(clusters.get(index).getData() > stCluster.getData()){
						clusters.get(index).setData(stCluster.getData());
					}
				}
			}
		}

		//order the clusters by size decreasing
		for(int i=1;i<clusters.size();i++){
			STITreeCluster<Double> cl1 = clusters.get(i);
			for(int j=0;j<i;j++){
				STITreeCluster<Double> cl2 = clusters.get(j);
				if(cl2.getClusterSize() < cl1.getClusterSize()){
					clusters.remove(cl1);
					clusters.add(j,cl1);
					break;
				}
			}
		}

		//adjust the time
		for(int i=0;i<clusters.size();i++){
			STITreeCluster<Double> cl1 = clusters.get(i);
			for(int j=i+1;j<clusters.size();j++){
				STITreeCluster<Double> cl2 = clusters.get(j);
				if(cl1.containsCluster(cl2)){
					if(cl1.getData() < cl2.getData()){
						cl2.setData(cl1.getData());
					}
				}
			}
		}

		//Add clusters of size 1.
		for (String t : stTaxa) {
			STITreeCluster<Double> c = new STITreeCluster<Double>(stTaxa);
			c.addLeaf(t);
			c.setData(0.0);
			clusters.add(c);
		}


		return clusters;
	}


	/**
	 * Computes the height of every node in a given tree
	 *
	 * @param 	tr	a given tree
	 */
	private void setNodeData(Tree tr){
		for(TNode n : tr.postTraverse()){
			if(n.isLeaf()){
				((STINode<Double>)n).setData(0.0);
			}
			else{
				for(TNode n_ch : n.getChildren()){
					double time = ((STINode<Double>)n_ch).getData()+n_ch.getParentDistance();
					if(time==TNode.NO_DISTANCE){
						throw new RuntimeException("Gene trees have branches that don't have time");
					}
					((STINode<Double>)n).setData(time);
					break;
				}
			}
		}
	}

	class Branch{
		STITreeCluster<Double> _child;
		STITreeCluster<Double> _parent;
		double _el;
		double _min_cost;
		double _max_score;
		Branch _min_lb;			// Left cluster for the min tree.
		Branch _min_rb;			// Right cluster for the min tree.
		List<Branch> _subBrs;


		public Branch(STITreeCluster<Double> child,STITreeCluster<Double> parent){
			_child = child;
			_parent = parent;
			_el = 0;
			_min_cost = -1;
			_max_score = -1;
			_min_lb = null;
			_min_rb = null;
		}

		public STITreeCluster<Double> getChildCluster(){
			return _child;
		}

		public STITreeCluster<Double> getParentCluster(){
			return _parent;
		}

		public BitSet getParentBitSet(){
			return _parent.getCluster();
		}

		public void setExtraLineage(double num){
			_el = num;
		}

		public void setMinCost(double cost){
			_min_cost = cost;
		}

		public double getExtraLineage(){
			return _el;
		}

		public double getMinCost(){
			return _min_cost;
		}

		public void setMaxScore(double score){
			_max_score = score;
		}

		public double getMaxScore(){
			return _max_score;
		}

		public String[] getTaxa(){
			return _child.getTaxa();
		}

		public void setMinLeftBranch(Branch br){
			_min_lb = br;
		}

		public void setMinRightBranch(Branch br){
			_min_rb = br;
		}

		public Branch getMinRightBranch(){
			return _min_rb;
		}

		public Branch getMinLeftBranch(){
			return _min_lb;
		}

		public boolean containsBranch(Branch br){
			if(this._child.containsCluster(br.getParentCluster())){
				return true;
			}
			else{
				return false;
			}
		}

		public boolean isCompatible(Branch b){
			if(this._parent.isDisjoint(b.getParentCluster())){
				return true;
			}
			else if(this._parent.equals(b.getParentCluster())){
				if(this._child.isDisjoint(b.getChildCluster())){
					return true;
				}
				if(this._child.equals(b.getChildCluster())){
					return false;
				}
			}
			else if((this._parent.containsCluster(b.getParentCluster())) &&
					(this._parent.getData() >= b.getParentCluster().getData())){
				if((this._child.containsCluster(b.getParentCluster())) &&
						(this._child.getData() >= b.getParentCluster().getData())){
					return true;
				}
				else if(this._child.isDisjoint(b.getParentCluster())){
					return true;
				}
			}
			else if((b.getParentCluster().containsCluster(this._parent)) &&
					(this._parent.getData() <= b.getParentCluster().getData())){
				if((b.getChildCluster().containsCluster(this._parent)) &&
						(this._parent.getData() <= b.getChildCluster().getData())){
					return true;
				}
				else if(this._parent.isDisjoint(b.getChildCluster())){
					return true;
				}
			}
			return false;
		}

		public String toString(){
			StringBuffer out = new StringBuffer();
			out.append(_child.toString());
			out.append("~");
			if(_parent!=null){
				out.append(_parent.toString());
			}
			out.append("/");
			out.append("("+_el+","+_min_cost+")");

			return out.toString();
		}

		public boolean equals(Object o){
			Branch b = (Branch)o;
			if(b.getParentCluster().equals(_parent) && b.getChildCluster().equals(_child)){
				return true;
			}
			else{
				return false;
			}
		}
	}

}
