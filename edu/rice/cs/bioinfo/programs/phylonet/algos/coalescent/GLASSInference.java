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

import edu.rice.cs.bioinfo.programs.phylonet.structs.sequence.model.SequenceAlignment;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/14/11
 * Time: 1:40 PM
 * To change this template use File | Settings | File Templates.
 */
public class GLASSInference
{
    private static boolean _print = true;

	private enum InputType {
		TREE,
		MATRIX
	}


	public static TaxaDistanceMatrix getTaxaDistanceMatrix(File inputFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(inputFile));
		String line = br.readLine().trim();
		int taxon_num = Integer.parseInt(line);
		line = br.readLine();
		String[] taxa = line.trim().split(" ");
		TaxaDistanceMatrix distanceMatrix = new TaxaDistanceMatrix(taxa);
		if (taxa.length != taxon_num) {
			throw new RuntimeException("The input file is in wrong format");
		}
		for (int i=0; i<taxon_num-1; i++) {
			line = br.readLine();
			String[] dis = line.trim().split(" ");
			if (dis.length != (taxon_num-1-i)) {
				throw new RuntimeException("The input file is in wrong format");
			}
			int pos = 0;
			for (int j=i+1; j<taxon_num; j++) {
				STITreeCluster c = new STITreeCluster(taxa);
				c.addLeaf(taxa[i]);
				c.addLeaf(taxa[j]);
				distanceMatrix.put(c, Double.parseDouble(dis[pos++]));
			}
		}
		while ((line=br.readLine())!=null) {
			if (line.trim().length() != 0) {
				throw new RuntimeException("The input file is in wrong format");
			}
		}
		return distanceMatrix;
		// for sequences
		/*
		br = new BufferedReader(new FileReader(inputFile));
		sequences = new ArrayList<SequenceAlignment>();
		try{
			while ((line = br.readLine()) != null) {
				if(line.trim().equals("start;") || line.trim().equals("start")){
					List<String> fileContent = new ArrayList<String>();
					while(true){
						line = br.readLine().trim();
						if(line == null){
							throw new SequenceException("Keyword end expected");
						}
						if(line.equals("end")||line.equals("end;")){
							break;
						}
						fileContent.add(line);
					}
					SequenceAlignment sa = new SequenceAlignment();
					sa.readPhyloNetSequences(fileContent);
					sequences.add(sa);
				}
			}
		}catch(SequenceException e){
			System.err.println(e);
			e.getStackTrace();
		}
		if (sequences.size() == 0) {
			System.err.println("Empty list of sequences. This may be caused by wrong input file format. " );
			System.err.println("The function exits.");
			return;
		}
		*/
	}



	/**
	 * Infers the species tree from the given list of gene trees with multiple alleles
	 *
	 * @param 	trees	the list of gene trees from which the species tree is to be
	 * 					inferred
	 * @param	taxonMap	Maps gene tree taxa to species tree taxa.
	 * @return	the tree inferred
	 *
	 */
	public Tree inferSpeciesTree(List<Tree> trees, Map<String,String> taxonMap){
		String error = Trees.checkMapping(trees, taxonMap);
		if(error!=null){
			throw new RuntimeException("Gene trees have leaf named " + error + " that hasn't been defined in the mapping file");
		}


		LinkedList<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		ArrayList<STITreeCluster> allClusters = new ArrayList<STITreeCluster>();

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

		for(Tree tr:trees){
			setNodeTime(tr);
		}

		Map<STITreeCluster,Double> taxonPairTime = computeTaxonPairTime(trees,stTaxa,gtTaxa,taxonMap);

		//System.out.println(taxonPairTime);
		for (String taxon: stTaxa){
			STITreeCluster cluster = new STITreeCluster(stTaxa);
			cluster.addLeaf(taxon);
			allClusters.add(cluster);
		}
		while(allClusters.size()>2)
			clutering(minClusters,allClusters,stTaxa,taxonPairTime);

		Tree inferredTree = Trees.buildTreeFromClusters(minClusters);
		return inferredTree;
	}


	/**
	 * Infers the species tree from the given list of gene trees with single alleles
	 *
	 * @param 	trees	the list of gene trees from which the species tree is to be inferred
	 * @return	the tree inferred
	 *
	 */
	public Tree inferSpeciesTree(List<Tree> trees){
		LinkedList<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		ArrayList<STITreeCluster> allClusters = new ArrayList<STITreeCluster>();

		List<String> temp = new LinkedList<String>();
		Tree t1 = null;
		for(Tree tr:trees){
			if(t1 == null){
				t1 = tr;
				for(TNode node:tr.getRoot().getLeaves()){
					temp.add(node.getName());
				}
			}
			else{
				if(!Trees.leafSetsAgree(t1, tr)){
					for(TNode node:tr.getRoot().getLeaves()){
						if(!temp.contains(node.getName())){
								temp.add(node.getName());
						}
					}
				}
			}
		}

		String[] taxa = new String[temp.size()];
		int index = 0;
		for(String taxon:temp){
			taxa[index++] = taxon;
		}

		for(Tree tr:trees){
			setNodeTime(tr);
		}

		Map<STITreeCluster,Double> taxonPairTime = computeTaxonPairTime(trees,taxa);

		//System.out.println(taxonPairTime);
		for (String taxon: taxa){
			STITreeCluster cluster = new STITreeCluster(taxa);
			cluster.addLeaf(taxon);
			allClusters.add(cluster);
		}

		while(allClusters.size()>2)
			clutering(minClusters,allClusters,taxa,taxonPairTime);

		Tree inferredTree = Trees.buildTreeFromClusters(minClusters);
		return inferredTree;
	}


	/**
	 * Infers the species tree from the given distances of taxon pairs
	 * @param 	distanceMatrix	taxa in species tree
	 * @return	the tree inferred
	 */
	public Tree inferSpeciesTreeFromTaxa(TaxaDistanceMatrix distanceMatrix) {
		String[] taxa = distanceMatrix.getTaxa();
		Map<STITreeCluster,Double> taxonPairTime = distanceMatrix.getDistanceMap();
		LinkedList<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		ArrayList<STITreeCluster> allClusters = new ArrayList<STITreeCluster>();

		for (String taxon: taxa){
			STITreeCluster cluster = new STITreeCluster(taxa);
			cluster.addLeaf(taxon);
			allClusters.add(cluster);
		}

		while(allClusters.size()>2)
			clutering(minClusters,allClusters,taxa,taxonPairTime);

		Tree inferredTree = Trees.buildTreeFromClusters(minClusters);
		return inferredTree;
	}


	public Tree inferSpeciesTree2(List<SequenceAlignment> sequences){
		LinkedList<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		ArrayList<STITreeCluster> allClusters = new ArrayList<STITreeCluster>();

		List<String> taxaList = new LinkedList<String>();

		for(SequenceAlignment sa: sequences){
			String[] taxa = sa.getTaxa();
			for(String taxon: taxa){
				if(!taxaList.contains(taxon)){
					taxaList.add(taxon);
				}
			}
		}

		String[] taxa = new String[taxaList.size()];
		for(int i = 0; i < taxa.length; i++){
			taxa[i] = taxaList.get(i);
		}

		Map<STITreeCluster,Double> clusterPairDistance = computeSequencePairDistance(sequences,taxa);

		for (String taxon: taxa){
			STITreeCluster cluster = new STITreeCluster(taxa);
			cluster.addLeaf(taxon);
			allClusters.add(cluster);
		}

		while(allClusters.size()>2)
			clutering(minClusters,allClusters,taxa,clusterPairDistance);

		Tree inferredTree = Trees.buildTreeFromClusters(minClusters);
		return inferredTree;
	}

	public Tree inferSpeciesTree2(List<SequenceAlignment> sequences,double percentage){
		LinkedList<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		ArrayList<STITreeCluster> allClusters = new ArrayList<STITreeCluster>();

		List<String> taxaList = new LinkedList<String>();

		for(SequenceAlignment sa: sequences){
			String[] taxa = sa.getTaxa();
			for(String taxon: taxa){
				if(!taxaList.contains(taxon)){
					taxaList.add(taxon);
				}
			}
		}

		String[] taxa = new String[taxaList.size()];
		for(int i = 0; i < taxa.length; i++){
			taxa[i] = taxaList.get(i);
		}

		Map<STITreeCluster,Double> clusterPairDistance = computeSequencePairDistance(sequences,taxa,percentage);

		for (String taxon: taxa){
			STITreeCluster cluster = new STITreeCluster(taxa);
			cluster.addLeaf(taxon);
			allClusters.add(cluster);
		}

		while(allClusters.size()>2)
			clutering(minClusters,allClusters,taxa,clusterPairDistance);

		Tree inferredTree = Trees.buildTreeFromClusters(minClusters);
		return inferredTree;
	}

	//TODO
	public Tree inferSpeciesTree2(List<SequenceAlignment> sequences, Map<String,String> taxonMap){
		LinkedList<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		ArrayList<STITreeCluster> allClusters = new ArrayList<STITreeCluster>();

		List<String> temp = new ArrayList<String>();
		for (String s : taxonMap.keySet()) {
			if (!temp.contains(taxonMap.get(s))) {
				temp.add(taxonMap.get(s));	// Species tree taxa.
			}
		}

		String[] stTaxa = new String[temp.size()];
		for (int i = 0; i < stTaxa.length; i++) {
			stTaxa[i] = temp.get(i);
		}


		Map<STITreeCluster,Double> clusterPairDistance = computeSequencePairDistance(sequences,stTaxa,taxonMap);


		//System.out.println(clusterPairDistance);
		for (String taxon: stTaxa){
			STITreeCluster cluster = new STITreeCluster(stTaxa);
			cluster.addLeaf(taxon);
			allClusters.add(cluster);
		}

		while(allClusters.size()>2)
			clutering(minClusters,allClusters,stTaxa,clusterPairDistance);

		Tree inferredTree = Trees.buildTreeFromClusters(minClusters);
		return inferredTree;
	}

	public Tree inferSpeciesTree2(List<SequenceAlignment> sequences, Map<String,String> taxonMap,double percentage){
		LinkedList<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		ArrayList<STITreeCluster> allClusters = new ArrayList<STITreeCluster>();

		List<String> temp = new ArrayList<String>();
		for (String s : taxonMap.keySet()) {
			if (!temp.contains(taxonMap.get(s))) {
				temp.add(taxonMap.get(s));	// Species tree taxa.
			}
		}

		String[] stTaxa = new String[temp.size()];
		for (int i = 0; i < stTaxa.length; i++) {
			stTaxa[i] = temp.get(i);
		}


		Map<STITreeCluster,Double> clusterPairDistance = computeSequencePairDistance(sequences,stTaxa,taxonMap,percentage);


		//System.out.println(clusterPairDistance);
		for (String taxon: stTaxa){
			STITreeCluster cluster = new STITreeCluster(stTaxa);
			cluster.addLeaf(taxon);
			allClusters.add(cluster);
		}

		while(allClusters.size()>2)
			clutering(minClusters,allClusters,stTaxa,clusterPairDistance);

		Tree inferredTree = Trees.buildTreeFromClusters(minClusters);
		return inferredTree;
	}


	private double computeHammingDistance(String seq1, String seq2){
		double distance = 0;
		//double B = 1-(0.3*0.3+0.2*0.2+0.2*0.2+0.3*0.3);
		double B = 0.75;
		char[] seq1_array = seq1.toCharArray();
		char[] seq2_array = seq2.toCharArray();
		for(int i=0;i<seq1_array.length;i++){
			if(!(seq1_array[i] == seq2_array[i])){
				distance++;
			}
		}
		//return distance;
		distance = distance / seq1_array.length;
		distance = (-1)*B*Math.log(1-distance/B);
		return distance;
	}



	/**
	 * This function is written to compute the distance of taxon pairs from a set of gene trees with multiple alleles
	 *
	 * @param 	trees	the list of gene trees
	 * @param	stTaxa	taxa in species tree
	 * @param	gtTaxa	taxa in gene trees
	 * @param	taxonMap	Maps gene tree taxa to species tree taxa.
	 * @return  maps from taxon pairs to distance
	 */
	private Map<STITreeCluster,Double> computeTaxonPairTime(List<Tree> trees, String[] stTaxa, String[] gtTaxa, Map<String,String> taxonMap){
		LinkedList<STITreeCluster<Double>> clusters = new LinkedList<STITreeCluster<Double>>();

		//Add the cluster containing all taxa to the end of the list.
		STITreeCluster<Double> all = new STITreeCluster<Double>(stTaxa);
		for (String t : stTaxa) {
			all.addLeaf(t);
		}
		double minTime = -1;
		for(Tree tr:trees){
			TNode root = tr.getRoot();
			if(((STINode<Double>)root).getData() < minTime || minTime == -1){
				minTime = ((STINode<Double>)root).getData();
			}
		}
		all.setData(minTime);

		clusters.add(all);


		for (Tree tr : trees) {
			for (STITreeCluster<Double> tc : tr.getClusters(gtTaxa, true)) {
				STITreeCluster<Double> stCluster = new STITreeCluster<Double>(stTaxa);
				for (String s : tc.getClusterLeaves()) {
					stCluster.addLeaf(taxonMap.get(s));
				}
				stCluster.setData(tc.getData());

				if(!clusters.contains(stCluster)){
					if(stCluster.getClusterSize()>1){
						clusters.add(stCluster);
					}
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


		Map<STITreeCluster,Double> taxonPairTime = new HashMap<STITreeCluster, Double>();
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

		//System.out.println(clusters);

		for(STITreeCluster<Double> cl1: clusters){
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
	 * This function is written to compute the distance of taxon pairs from a set of gene trees with single allele
	 *
	 * @param 	trees	the list of gene trees
	 * @param	taxa	taxa in species tree/gene trees
	 * @return  maps from taxon pairs to distance
	 */
	private Map<STITreeCluster,Double> computeTaxonPairTime(List<Tree> trees, String[] taxa){
		LinkedList<STITreeCluster<Double>> clusters = new LinkedList<STITreeCluster<Double>>();

		//Add the cluster containing all taxa to the first of the list.
		STITreeCluster<Double> all = new STITreeCluster<Double>(taxa);
		for (String t : taxa) {
			all.addLeaf(t);
		}
		double minTime = -1;
		for(Tree tr:trees){
			TNode root = tr.getRoot();
			if(((STINode<Double>)root).getData() < minTime || minTime == -1){
				minTime = ((STINode<Double>)root).getData();
			}
		}
		all.setData(minTime);

		clusters.add(all);


		for (Tree tr : trees) {
			for (STITreeCluster<Double> tc : tr.getClusters(taxa, true)) {
				STITreeCluster<Double> stCluster = new STITreeCluster<Double>(tc);
				stCluster.setData(tc.getData());

				if(!clusters.contains(stCluster)){
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


		Map<STITreeCluster,Double> taxonPairTime = new HashMap<STITreeCluster, Double>();
		ArrayList<STITreeCluster<Double>> taxonPairs = new ArrayList<STITreeCluster<Double>>();

		for(int i=0;i<taxa.length;i++){
			for(int j=i+1;j<taxa.length;j++){
				STITreeCluster<Double> cl = new STITreeCluster<Double>(taxa);
				BitSet bs = new BitSet(taxa.length);
				bs.set(i);
				bs.set(j);
				cl.setCluster(bs);
				cl.setData(-1.0);
				taxonPairs.add(cl);
			}
		}

		//System.out.println(clusters);

		for(STITreeCluster<Double> cl1: clusters){
			for(STITreeCluster<Double> cl2 : taxonPairs){
				if(cl1.containsCluster(cl2)){
					if(cl1.getData()<cl2.getData() || cl2.getData()==-1){
						cl2.setData(cl1.getData());
					}
				}
			}
		}

		for(STITreeCluster<Double> clwt : taxonPairs){
			STITreeCluster cl = new STITreeCluster(taxa);
			cl.setCluster(clwt.getCluster());
			taxonPairTime.put(cl, clwt.getData());
		}

		return taxonPairTime;
	}

	private Map<STITreeCluster,Double> computeSequencePairDistance(List<SequenceAlignment> sequences, String[] stTaxa){
		Map<STITreeCluster,Double> sequencePairDistance = new HashMap<STITreeCluster,Double>();

		for(SequenceAlignment sa : sequences){
			String[] taxa = sa.getTaxa();
			for(int i=0;i<taxa.length;i++){
				String seq1 = sa.getSequences()[i];
				for(int j=i+1;j<taxa.length;j++){
					String seq2 = sa.getSequences()[j];
					Double distance = computeHammingDistance(seq1,seq2);
					STITreeCluster cluster = new STITreeCluster(stTaxa);
					cluster.addLeaf(taxa[i]);
					cluster.addLeaf(taxa[j]);
					Double old_distance = sequencePairDistance.get(cluster);
					if(old_distance == null || old_distance > distance){
						sequencePairDistance.put(cluster, distance);
					}
				}
			}
		}
		return sequencePairDistance;
	}

	//for remove bottom sequences
	private Map<STITreeCluster,Double> computeSequencePairDistance(List<SequenceAlignment> sequences, String[] stTaxa,double percentage){
		Map<STITreeCluster,List<Double>> sequencePairDistances = new HashMap<STITreeCluster,List<Double>>();

		for(SequenceAlignment sa : sequences){
			String[] taxa = sa.getTaxa();
			for(int i=0;i<taxa.length;i++){
				String seq1 = sa.getSequences()[i];
				for(int j=i+1;j<taxa.length;j++){
					String seq2 = sa.getSequences()[j];
					Double distance = computeHammingDistance(seq1,seq2);
					STITreeCluster cluster = new STITreeCluster(stTaxa);
					cluster.addLeaf(taxa[i]);
					cluster.addLeaf(taxa[j]);
					List<Double> dist_list = sequencePairDistances.get(cluster);
					if(dist_list == null){
						List<Double> newlist = new ArrayList<Double>();
						newlist.add(distance);
						sequencePairDistances.put(cluster, newlist);
					}
					else{
						boolean find = false;
						for(int k=0;k<dist_list.size();k++){
							if(dist_list.get(k)>distance){
								dist_list.add(k, distance);
								find = true;
								break;
							}
						}
						if(!find){
							dist_list.add(distance);
						}
					}
				}
			}
		}


		Map<STITreeCluster,Double> sequencePairDistance = new HashMap<STITreeCluster,Double>();
		for(Map.Entry<STITreeCluster, List<Double>> entry: sequencePairDistances.entrySet()){
			List<Double> dist_list = entry.getValue();
			int eliminate = dist_list.size()-(int)(dist_list.size()*(1-percentage));
			if(eliminate>=dist_list.size()){
				eliminate=dist_list.size()-1;
			}
			sequencePairDistance.put(entry.getKey(), dist_list.get(eliminate));
		}

		return sequencePairDistance;
	}

	private Map<STITreeCluster,Double> computeSequencePairDistance(List<SequenceAlignment> sequences,String[] stTaxa,Map<String,String> taxonMap){
		Map<STITreeCluster,Double> sequencePairDistance = new HashMap<STITreeCluster,Double>();
		for(SequenceAlignment sa : sequences){
			String[] gtTaxa = sa.getTaxa();
			for(int i=0;i<gtTaxa.length;i++){
				String seq1 = sa.getSequences()[i];
				String s1 = taxonMap.get(gtTaxa[i]);
				for(int j=i+1;j<gtTaxa.length;j++){
					String s2 = taxonMap.get(gtTaxa[j]);
					if(s1.equals(s2)){
						continue;
					}
					String seq2 = sa.getSequences()[j];
					Double distance = computeHammingDistance(seq1,seq2);
					STITreeCluster cluster = new STITreeCluster(stTaxa);
					cluster.addLeaf(s1);
					cluster.addLeaf(s2);
					Double old_distance = sequencePairDistance.get(cluster);
					if(old_distance == null || old_distance > distance){
						sequencePairDistance.put(cluster, distance);
					}
				}
			}
		}
		return sequencePairDistance;
	}

	private Map<STITreeCluster,Double> computeSequencePairDistance(List<SequenceAlignment> sequences,String[] stTaxa,Map<String,String> taxonMap,double percentage){
		Map<STITreeCluster,List<Double>> sequencePairDistances = new HashMap<STITreeCluster,List<Double>>();
		for(SequenceAlignment sa : sequences){
			String[] gtTaxa = sa.getTaxa();
			for(int i=0;i<gtTaxa.length;i++){
				String seq1 = sa.getSequences()[i];
				String s1 = taxonMap.get(gtTaxa[i]);
				for(int j=i+1;j<gtTaxa.length;j++){
					String s2 = taxonMap.get(gtTaxa[j]);
					if(s1.equals(s2)){
						continue;
					}
					String seq2 = sa.getSequences()[j];
					Double distance = computeHammingDistance(seq1,seq2);
					STITreeCluster cluster = new STITreeCluster(stTaxa);
					cluster.addLeaf(s1);
					cluster.addLeaf(s2);
					List<Double> dist_list = sequencePairDistances.get(cluster);
					if(dist_list == null){
						List<Double> newlist = new ArrayList<Double>();
						newlist.add(distance);
						sequencePairDistances.put(cluster, newlist);
					}
					else{
						boolean find = false;
						for(int k=0;k<dist_list.size();k++){
							if(dist_list.get(k)>distance){
								dist_list.add(k, distance);
								find = true;
								break;
							}
						}
						if(!find){
							dist_list.add(distance);
						}
					}
				}
			}
		}

		Map<STITreeCluster,Double> sequencePairDistance = new HashMap<STITreeCluster,Double>();
		for(Map.Entry<STITreeCluster, List<Double>> entry: sequencePairDistances.entrySet()){
			List<Double> dist_list = entry.getValue();
			int eliminate = dist_list.size()-(int)(dist_list.size()*(1-percentage));
			if(eliminate>=dist_list.size()){
				eliminate=dist_list.size()-1;
			}
			sequencePairDistance.put(entry.getKey(), dist_list.get(eliminate));
		}

		return sequencePairDistance;
	}


	private Map<STITreeCluster,Double> computeSequencePairDistance2(List<SequenceAlignment> sequences,String[] stTaxa,Map<String,String> taxonMap,double percentage){
		ArrayList<STITreeCluster<Double>> allPairDistance = new ArrayList<STITreeCluster<Double>>();
		for(SequenceAlignment sa : sequences){
			String[] gtTaxa = sa.getTaxa();
			for(int i=0;i<gtTaxa.length;i++){
				String seq1 = sa.getSequences()[i];
				String s1 = taxonMap.get(gtTaxa[i]);
				for(int j=i+1;j<gtTaxa.length;j++){
					String s2 = taxonMap.get(gtTaxa[j]);
					if(s1.equals(s2)){
						continue;
					}
					String seq2 = sa.getSequences()[j];
					Double distance = computeHammingDistance(seq1,seq2);
					STITreeCluster<Double> cluster = new STITreeCluster<Double>(stTaxa);
					cluster.addLeaf(s1);
					cluster.addLeaf(s2);
					cluster.setData(distance);
					allPairDistance.add(cluster);
				}
			}
		}
		return doEliminating(allPairDistance,stTaxa,percentage);
	}

	private Map<STITreeCluster,Double> doEliminating(ArrayList<STITreeCluster<Double>> allPairDistance,String[] stTaxa,double percentage){
		for(int i=1;i<allPairDistance.size();i++){
			STITreeCluster<Double> clwd1 = allPairDistance.get(i);
			for(int j=0;j<i;j++){
				STITreeCluster<Double> clwd2 = allPairDistance.get(j);
				if(clwd1.getData() < clwd2.getData()){
					allPairDistance.remove(i);
					allPairDistance.add(j,clwd1);
					break;
				}
			}
		}

		if(allPairDistance.size()>(stTaxa.length*(stTaxa.length-1)/2)){
			int eliminate = allPairDistance.size()-(int)(allPairDistance.size()*(1-percentage));
			for(int i=0;i<eliminate;i++){
				allPairDistance.remove(0);
			}
		}


		Map<STITreeCluster,Double> sequencePairDistance = new HashMap<STITreeCluster,Double>();
		for(int i=0;i<stTaxa.length;i++){
			for(int j=i+1;j<stTaxa.length;j++){
				STITreeCluster cl = new STITreeCluster(stTaxa);
				cl.addLeaf(stTaxa[i]);
				cl.addLeaf(stTaxa[j]);
				int index = allPairDistance.indexOf(cl);
				Double distance;
				if(index==-1){
					distance = -1.0;
				}
				else{
					distance = allPairDistance.get(index).getData();
				}
				sequencePairDistance.put(cl, distance);
			}
		}
		return sequencePairDistance;
	}



	/**
	 * This function is written to compute the depth of every node in a given tree
	 *
	 * @param 	tr	a given tree
	 */
	private void setNodeTime(Tree tr){
		for(TNode n : tr.postTraverse()){
			if(n.isLeaf()){
				((STINode<Double>)n).setData(0.0);
			}
			else{
				for(TNode n_ch : n.getChildren()){
					Double time = ((STINode<Double>)n_ch).getData()+n_ch.getParentDistance();
					((STINode<Double>)n).setData(time);
					break;
				}
			}
		}
	}


	/**
	 * This function is written to merge two clusters having the smallest distance
	 *
	 * @param 	minClusters		a list of the clusters merged
	 * @param 	allClusters		a list of all clusters using during merging
	 * @param 	stTaxa		taxa in species tree
	 * @param 	taxonPairTime		maps from taxon pair to distance
	 */
	private void clutering(List<STITreeCluster> minClusters,List<STITreeCluster> allClusters,String[] stTaxa,Map<STITreeCluster,Double> taxonPairTime){
		double minTime=-1;
		int mergingClusterA=0, mergingClusterB=0;

		for(int i=0; i<allClusters.size(); i++){
			STITreeCluster clusterA = allClusters.get(i);
			for(int j=i+1; j<allClusters.size(); j++){
				STITreeCluster clusterB = allClusters.get(j);
				for (String taxonA: clusterA.getClusterLeaves())
					for (String taxonB: clusterB.getClusterLeaves()){
						STITreeCluster cl = new STITreeCluster(stTaxa);
						cl.addLeaf(taxonA);
						cl.addLeaf(taxonB);
						double time = taxonPairTime.get(cl);
						if(minTime==-1 || minTime>time){
							minTime = time;
							mergingClusterA = i;
							mergingClusterB = j;
						}
					}
			}
		}

		allClusters.get(mergingClusterA).getCluster().or(allClusters.get(mergingClusterB).getCluster());
		allClusters.remove(mergingClusterB);
		STITreeCluster newCluster = new STITreeCluster(stTaxa);
		BitSet b = new BitSet(stTaxa.length);
		b.or(allClusters.get(mergingClusterA).getCluster());
		newCluster.setCluster(b);
		minClusters.add(newCluster);
	}


}
