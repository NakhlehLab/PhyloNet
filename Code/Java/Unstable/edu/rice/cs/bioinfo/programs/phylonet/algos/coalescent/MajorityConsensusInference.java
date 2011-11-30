package edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeClusterWD;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/15/11
 * Time: 1:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class MajorityConsensusInference
{

	/**
	 * Returns a taxon map constructed from the contents of the specified file.
	 *
	 * @param file			the file whose contents are to be parsed
	 * @return				a String-to-String allele mapping
	 * @throws IOException 	if the supplied file could not be opened or read
	 */
	public static Map<String, String> getTaxonMapFromFile(File file) throws IOException {
		Map<String, String> taxonMap = new HashMap<String,String>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		try {
			String line;
			while ((line = br.readLine()) != null) {
				String[] mapString = line.trim().split(";");
				for(String s : mapString) {
					String species = s.substring(0,s.indexOf(":")).trim();
					s = s.substring(s.indexOf(":") + 1);
					String[] alleles = s.split(",");
					for(String allele:alleles) {
						allele = allele.trim();
						if (taxonMap.containsKey(allele))
							throw new RuntimeException("an allele can only map to one species");
						taxonMap.put(allele, species);
					}
				}
			}
		} finally {
			br.close();
		}
		return taxonMap;
	}

	public Tree inferSpeciesTree(List<Tree> trees, boolean rooted){
		if(rooted){
			return inferSpeciesTreeRooted(trees);
		}
		else{
			return inferSpeciesTreeUnrooted(trees);
		}
	}

	public Tree inferSpeciesTree(List<Tree> trees, boolean rooted, Map<String, String> taxonMap){
		String error = Trees.checkMapping(trees, taxonMap);
		if(error!=null){
			throw new RuntimeException("Gene trees have leaf named " + error + "that hasn't been defined in the mapping file");
		}

		if(rooted){
			return inferSpeciesTreeRooted(trees, taxonMap);
		}
		else{
			return inferSpeciesTreeUnrooted(trees, taxonMap);
		}
	}

	/**
	 * Infers the species tree from the given list of rooted gene trees with single allele.
	 *
	 * @param 	trees	the list of gene trees from which the species tree is to be inferred
	 *
	 * @return	inferred species tree
	 */
    public Tree inferSpeciesTreeRooted(List<Tree> trees){

		// make sure that all binary nodes are removed (to avoid double counting edges)
		// make sure that trees have the same leaf set

		for(Tree tr : trees) {
			Trees.removeBinaryNodes(new STITree<Object>(tr));
		}


		List<String> taxalist = new ArrayList<String>();
		for(Tree tr: trees){
			for (TNode node : tr.postTraverse()) {
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


		LinkedList<STITreeClusterWD<Integer>> cls = new LinkedList<STITreeClusterWD<Integer>>();

		for(Tree tr: trees){
			for(STITreeClusterWD<Integer> cl : tr.getClustersWD(taxa,true)){
				index = cls.indexOf(cl);
				if(index!=-1){
					STITreeClusterWD<Integer> clwd = cls.get(index);
					clwd.setData(clwd.getData()+1);
				}
				else{
					STITreeClusterWD<Integer> clwd = cl.duplicate();
					clwd.setData(1);
					cls.add(clwd);
				}
			}
		}

		//order the clusters by its frequency
		for(int i=1;i<cls.size();i++){
			STITreeClusterWD<Integer> cl1 = cls.get(i);
			for(int j=0;j<i;j++){
				STITreeClusterWD<Integer> cl2 = cls.get(j);
				if(cl2.getData() < cl1.getData()){
					cls.remove(cl1);
					cls.add(j,cl1);
					break;
				}
			}
		}

		//get the clusters for result trees
		List<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();

		int size = trees.size();
		for(STITreeClusterWD<Integer> cl: cls){
			if(cl.getData()/(double)size > 0.5){
				minClusters.add(cl);
			}
			else{
				boolean compatible = true;
				for(STITreeCluster c_ex : minClusters){
					if(!c_ex.isCompatible(cl)){
						compatible = false;
						break;
					}
				}
				if(compatible){
					minClusters.add(cl);
				}
			}
			if(minClusters.size() == taxa.length-2){
				break;
			}
		}

		Tree tr = Trees.buildTreeFromClusters(minClusters);
		return tr;
    }


	/**
	 * Infers the species tree from the given list of rooted gene trees with multiple alleles.
	 *
	 * @param 	trees	the list of gene trees from which the species tree is to be inferred
	 * @param	taxonMap	Maps gene tree taxa to species tree taxa.
	 *
	 * @return	inferred species tree
	 */
    public Tree inferSpeciesTreeRooted(List<Tree> trees, Map<String,String> taxonMap){

		// make sure that all binary nodes are removed (to avoid double counting edges)
		// make sure that trees have the same leaf set

		for(Tree tr : trees) {
			Trees.removeBinaryNodes(new STITree<Object>(tr));
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

		LinkedList<STITreeClusterWD<Double>> cls = new LinkedList<STITreeClusterWD<Double>>();

		for(Tree tr: trees){
			Map<STITreeCluster, LinkedList<STITreeCluster>> cmap = new HashMap<STITreeCluster, LinkedList<STITreeCluster>>();
			//get all clusters
			for(STITreeCluster cl : tr.getClusters(gtTaxa, true)){
				STITreeCluster stCluster = new STITreeCluster(stTaxa);
				for (String s : cl.getClusterLeaves()) {
					stCluster.addLeaf(taxonMap.get(s));
				}
				if(stCluster.getClusterSize()==1 || stCluster.getClusterSize()==stTaxa.length){
					continue;
				}
				if(cmap.containsKey(stCluster)){
					LinkedList<STITreeCluster> l = cmap.get(stCluster);
					int relation = 0;
					LinkedList<STITreeCluster> rl = new LinkedList<STITreeCluster>();
					for(STITreeCluster c : l){
						if(c.containsCluster(cl)){
							relation = -1;
							break;
						}
						else if(cl.containsCluster(c)){
							relation = 1;
							rl.add(c);
						}
					}
					if(relation!=-1){
						l.add(cl);
					}
					if(rl.size()>0){
						for(STITreeCluster rc:rl){
							l.remove(rc);
						}
					}
				}
				else{
					LinkedList<STITreeCluster> l = new LinkedList<STITreeCluster>();
					l.add(cl);
					cmap.put(stCluster, l);
				}
			}

			int[] indNum = new int[stTaxa.length];
			for(int i=0;i<indNum.length;i++){
				indNum[i] = 0;
			}
			for(TNode node : tr.getRoot().getLeaves()){
				indNum[temp2.indexOf(taxonMap.get(node.getName()))]++;
			}

			//compute cluster frequency inside every tree
			for(STITreeCluster stcl : cmap.keySet()){
				/*int totalNum = 0;
				for(String leaf : stcl.getClusterLeaves()){
					totalNum += indNum[temp2.indexOf(leaf)];
				}
				int appearNum = 0;
				LinkedList<STITreeCluster> gtcls = cmap.get(stcl);
				for(STITreeCluster c:gtcls){
					appearNum += c.getClusterSize();
				}
				double frequency = appearNum/(totalNum+0.0);*/
				int[] num = new int[stTaxa.length];
				for(int i=0;i<num.length;i++){
					num[i] = 0;
				}

				LinkedList<STITreeCluster> gtcls = cmap.get(stcl);
				for(STITreeCluster c:gtcls){
					String alleles[] = c.getClusterLeaves();
					for(String allele : alleles){
						num[temp2.indexOf(taxonMap.get(allele))]++;
					}
				}

				double frequency = 0.0;

				for(int i=0;i<indNum.length;i++){
					if(num[i]!=0){
						frequency += num[i]/(indNum[i]+0.0);
					}
				}
				frequency = frequency/stcl.getClusterSize();

				//if(frequency>1)System.out.println("wrong");


				if(cls.contains(stcl)){
					STITreeClusterWD<Double> cl_ex = cls.get(cls.indexOf(stcl));
					cl_ex.setData(cl_ex.getData()+frequency);
				}
				else{
					STITreeClusterWD<Double> newcl = new STITreeClusterWD<Double>(stcl);
					newcl.setData(frequency);
					cls.add(newcl);
				}
			}
		}



		//order the clusters by its frequency decreasing
		for(int i=1;i<cls.size();i++){
			STITreeClusterWD<Double> cl1 = cls.get(i);
			for(int j=0;j<i;j++){
				STITreeClusterWD<Double> cl2 = cls.get(j);
//				if((cl2.getData().equals(cl1.getData()) && cl2.getClusterSize() > cl1.getClusterSize())
				if(cl2.getData() < cl1.getData()){
					cls.remove(cl1);
					cls.add(j,cl1);
					break;
				}
			}
		}

		//System.out.println(cls);
		//get the clusters for result trees
		List<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();

		for(STITreeClusterWD<Double> cl: cls){
			boolean compatible = true;
			for(STITreeCluster c_ex : minClusters){
				if(!c_ex.isCompatible(cl)){
					compatible = false;
					break;
				}
			}
			if(compatible){
				minClusters.add(cl);
			}
			if(minClusters.size() == stTaxa.length-2){
				break;
			}
		}

		return Trees.buildTreeFromClusters(minClusters);
    }


	/**
	 * Infers the species tree from the given list of unrooted gene trees with single allele.
	 *
	 * @param 	trees	the list of gene trees from which the species tree is to be inferred
	 *
	 * @return	inferred species tree
	 */
    public Tree inferSpeciesTreeUnrooted(List<Tree> trees){
		if(trees.size()==1){
			return trees.get(0);
		}

		// make sure that all binary nodes are removed (to avoid double counting edges)
		// make sure that trees have the same leaf set

		for(Tree tr : trees) {
			Trees.removeBinaryNodes(new STITree<Object>(tr));
		}


		List<String> taxalist = new ArrayList<String>();
		for(Tree tr: trees){
			for (TNode node : tr.postTraverse()) {
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

		List<STITreeClusterWD<Integer>> clusters = new LinkedList<STITreeClusterWD<Integer>>();

		for(Tree tr: trees){
			for(STITreeCluster cl : tr.getBipartitionClusters(taxa, false)){
				if(cl.getClusterSize()<=1 || cl.getClusterSize()>taxa.length-1){
					continue;
				}
				int pos = clusters.indexOf(cl);
				if(pos == -1){
					STITreeClusterWD<Integer> nc = new STITreeClusterWD<Integer>(cl);
					nc.setData(1);
					clusters.add(nc);
				}
				else{
					STITreeClusterWD<Integer> nc = clusters.get(pos);
					nc.setData(nc.getData()+1);
				}
			}
		}


		//order the clusters by its frequency
		for(int i=1;i<clusters.size();i++){
			STITreeClusterWD<Integer> cl1 = clusters.get(i);
			for(int j=0;j<i;j++){
				STITreeClusterWD<Integer> cl2 = clusters.get(j);
				if(cl2.getData() < cl1.getData()){
					clusters.remove(cl1);
					clusters.add(j,cl1);
					break;
				}
			}
		}
	/*
		for(STITreeClusterWD<Integer> c: clusters){
			System.out.println(c);
		}
	 */
		//get the clusters for result trees
		List<STITreeCluster> minClusters = new ArrayList<STITreeCluster>();
		for(STITreeCluster c: clusters){
			boolean compatible = true;
			for(STITreeCluster ex_c: minClusters){
				if(!c.isCompatible(ex_c)){
						compatible = false;
						break;
					}
			}
			if(compatible){
				minClusters.add(c);
				/*
				if(minClusters.size()==taxa.length-2){
					break;
				}
				*/
			}
		}
		return Trees.buildTreeFromClusters(minClusters);
    }


	/**
	 * Infers the species tree from the given list of unrooted gene trees with multiple alleles.
	 *
	 * @param 	trees	the list of gene trees from which the species tree is to be inferred
	 * @param	taxonMap	Maps gene tree taxa to species tree taxa.
	 *
	 * @return	inferred species tree
	 */

    public Tree inferSpeciesTreeUnrooted(List<Tree> trees, Map<String,String> taxonMap){

		// make sure that all binary nodes are removed (to avoid double counting edges)
		// make sure that trees have the same leaf set
		for(Tree tr : trees) {
			Trees.removeBinaryNodes(new STITree<Object>(tr));
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

		List<STITreeClusterWD<Double>> clusters = new LinkedList<STITreeClusterWD<Double>>();

		for(Tree tr: trees){
			Map<STITreeCluster, List<STITreeCluster>> cmap = new HashMap<STITreeCluster, List<STITreeCluster>>();

			for(STITreeCluster gtcl: tr.getBipartitionClusters(gtTaxa, false)){
				STITreeCluster stcl = new STITreeCluster(stTaxa);
				for(String leaf: gtcl.getClusterLeaves()){
					stcl.addLeaf(taxonMap.get(leaf));
				}
				int stcl_size = stcl.getClusterSize();
				if(stcl_size>1 && stcl_size<stTaxa.length){
					if(cmap.containsKey(stcl)){
						cmap.get(stcl).add(gtcl);
					}
					else{
						List<STITreeCluster> l = new ArrayList<STITreeCluster>();
						l.add(gtcl);
						cmap.put(stcl, l);
					}
				}
			}

			int[] indNum = new int[stTaxa.length];
			for(int i=0;i<indNum.length;i++){
				indNum[i] = 0;
			}

			for(TNode node : tr.getRoot().getLeaves()){
				indNum[temp2.indexOf(taxonMap.get(node.getName()))]++;
			}

			//compute cluster frequency inside every tree
			for(STITreeCluster stcl : cmap.keySet()){
				int[] num = new int[stTaxa.length];
				for(int i=0;i<num.length;i++){
					num[i] = 0;
				}
				List<STITreeCluster> gtcls = cmap.get(stcl);
				STITreeCluster temp = null;
				for(STITreeCluster gtcl: gtcls){
					if(temp == null){
						temp = new STITreeCluster(gtcl);
					}
					else{
						temp = temp.merge(gtcl);
					}
				}
				for(String leaf: temp.getClusterLeaves()){
					String st_leaf = taxonMap.get(leaf);
					num[temp2.indexOf(st_leaf)]++;
				}
				double frequency = 0.0;
				for(int i=0;i<indNum.length;i++){
					if(num[i]!=0){
						frequency += num[i]/(indNum[i]+0.0);
					}
				}
				frequency = frequency/stcl.getClusterSize();
				int pos = clusters.indexOf(stcl);
				if(pos == -1){
					STITreeClusterWD<Double> cd = new STITreeClusterWD<Double>(stcl);
					cd.setData(frequency);
					clusters.add(cd);
				}
				else{
					STITreeClusterWD<Double> cd = clusters.get(pos);
					cd.setData(frequency + cd.getData());
				}
			}
		}

/*
		int[] changed = new int[clusters.size()];
		for(int i=0;i<changed.length;i++){
			changed[i] = 0;
		}
		int num_trees = trees.size();
		for(int i=0;i<clusters.size();i++){
			if(changed[i]==1)continue;
			STITreeClusterWD<Double> cd = clusters.get(i);
			changed[i] = 1;
			if(cd.getClusterSize() == stTaxa.length-1){
				cd.setData(cd.getData()+num_trees);
			}
			else{

				int pos = clusters.indexOf(cd.complementaryCluster());
				if(pos != -1){
					STITreeClusterWD<Double> cc = clusters.get(pos);
					cd.setData(cd.getData()+cc.getData());
					cc.setData(cd.getData());
					changed[pos] = 1;
				}
			}
		}
*/

		//order the clusters by its frequency
		for(int i=1;i<clusters.size();i++){
			STITreeClusterWD<Double> cl1 = clusters.get(i);
			for(int j=0;j<i;j++){
				STITreeClusterWD<Double> cl2 = clusters.get(j);
				if(cl2.getData() < cl1.getData()){
					clusters.remove(cl1);
					clusters.add(j,cl1);
					break;
				}
			}
		}
/*
		for(STITreeClusterWD<Double> c: clusters){
			System.out.println(c);
		}
*/
		//get the clusters for result trees
		List<STITreeCluster> minClusters = new ArrayList<STITreeCluster>();

		for(STITreeCluster stcl: clusters){
			boolean comp = true;
			for(STITreeCluster cl_ex:minClusters){
				if(!cl_ex.isCompatible(stcl)){
					comp = false;
					break;
				}
			}
			if(comp){
				minClusters.add(stcl);
			}
			if(minClusters.size() == stTaxa.length-2){
				break;
			}
		}
		/*
		for(STITreeCluster c: minClusters){
			System.out.println(c);
		}
		*/
		return Trees.buildTreeFromClusters(minClusters);
    }



}
