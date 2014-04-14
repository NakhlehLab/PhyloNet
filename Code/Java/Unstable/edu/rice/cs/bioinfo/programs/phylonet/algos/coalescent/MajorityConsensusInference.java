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
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
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

	public Tree inferSpeciesTree(List<MutableTuple<Tree,Double>> trees, boolean rooted, int percentage){
		if(rooted){
			return inferSpeciesTreeRooted(trees, percentage);
		}
		else{
			return inferSpeciesTreeUnrooted(trees, percentage);
		}
	}

	public Tree inferSpeciesTree(List<MutableTuple<Tree,Double>> trees, boolean rooted, Map<String, String> taxonMap, int percentage){
		String error = Trees.checkMapping(trees, taxonMap);
		if(error!=null){
			throw new RuntimeException("Gene trees have leaf named " + error + " that hasn't been defined in the mapping file");
		}

		if(rooted){
			return inferSpeciesTreeRooted(trees, taxonMap, percentage);
		}
		else{
			return inferSpeciesTreeUnrooted(trees, taxonMap, percentage);
		}
	}

	/**
	 * Infers the species tree from the given list of rooted gene trees with single allele.
	 *
	 * @param 	trees	the list of gene trees from which the species tree is to be inferred
	 *
	 * @return	inferred species tree
	 */
    public Tree inferSpeciesTreeRooted(List<MutableTuple<Tree,Double>> trees, int percentage){
        double size = 0;
		Set<String> taxalist = new HashSet<String>();
		for(MutableTuple<Tree,Double> tr: trees){
			for (TNode node : tr.Item1.postTraverse()) {
				if (node.isLeaf()) {
					taxalist.add(node.getName());
				}
			}
            size += tr.Item2;
		}

		String[] taxa = new String[taxalist.size()];
		int index = 0;
		for(String taxon: taxalist){
			taxa[index++] = taxon;
		}


		LinkedList<STITreeCluster<Double>> cls = new LinkedList<STITreeCluster<Double>>();

        for(MutableTuple<Tree,Double> tuple: trees){
            Tree tr = tuple.Item1;
            double weight = tuple.Item2;
			for(STITreeCluster<Double> cl : tr.getClusters(taxa,false)){
				index = cls.indexOf(cl);
				if(index!=-1){
					STITreeCluster<Double> clwd = cls.get(index);
					clwd.setData(clwd.getData()+weight);
				}
				else{
					STITreeCluster<Double> clwd = cl.duplicate();
					clwd.setData(weight);
					cls.add(clwd);
				}
			}
		}

		//order the clusters by its frequency
		for(int i=1;i<cls.size();i++){
			STITreeCluster<Double> cl1 = cls.get(i);
			for(int j=0;j<i;j++){
				STITreeCluster<Double> cl2 = cls.get(j);
				if(cl2.getData() < cl1.getData()){
					cls.remove(cl1);
					cls.add(j,cl1);
					break;
				}
			}
		}

		//get the clusters for result trees
		List<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();

        if(percentage == 0){
            for(STITreeCluster<Double> cl: cls){
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
        }
        else{
            for(STITreeCluster<Double> cl: cls){
                if(cl.getData()/size > percentage/100.0){
                    minClusters.add(cl);
                }
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
    public Tree inferSpeciesTreeRooted(List<MutableTuple<Tree,Double>> trees, Map<String,String> taxonMap, int percentage){


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

		LinkedList<STITreeCluster<Double>> cls = new LinkedList<STITreeCluster<Double>>();

        double size = 0;
        for(MutableTuple<Tree,Double> tuple: trees){
            Tree tr = tuple.Item1;
            double weight = tuple.Item2;
            size += weight;
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
					STITreeCluster<Double> cl_ex = cls.get(cls.indexOf(stcl));
					cl_ex.setData(cl_ex.getData()+frequency);
				}
				else{
					STITreeCluster<Double> newcl = new STITreeCluster<Double>(stcl);
					newcl.setData(frequency);
					cls.add(newcl);
				}
			}
		}



		//order the clusters by its frequency decreasing
		for(int i=1;i<cls.size();i++){
			STITreeCluster<Double> cl1 = cls.get(i);
			for(int j=0;j<i;j++){
				STITreeCluster<Double> cl2 = cls.get(j);
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

		if(percentage == 0){
            for(STITreeCluster<Double> cl: cls){
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
        }

        else{
            for(STITreeCluster<Double> cl: cls){
                if(cl.getData()/size > percentage/100.0){
                    minClusters.add(cl);
                }
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
    public Tree inferSpeciesTreeUnrooted(List<MutableTuple<Tree,Double>> trees, int percentage){

        double size = 0;
        Set<String> taxalist = new HashSet<String>();
        for(MutableTuple<Tree,Double> tr: trees){
            for (TNode node : tr.Item1.postTraverse()) {
                if (node.isLeaf()) {
                    taxalist.add(node.getName());
                }
            }
            size += tr.Item2;
        }

		String[] taxa = new String[taxalist.size()];
		int index = 0;
		for(String taxon: taxalist){
			taxa[index++] = taxon;
		}

		List<STITreeCluster<Double>> clusters = new LinkedList<STITreeCluster<Double>>();

        for(MutableTuple<Tree,Double> tuple: trees){
            Tree tr = tuple.Item1;
            double weight = tuple.Item2;
			for(STITreeCluster cl : tr.getBipartitionClusters(taxa, false)){
				if(cl.getClusterSize()<=1 || cl.getClusterSize()>taxa.length-1){
					continue;
				}
				int pos = clusters.indexOf(cl);
				if(pos == -1){
					STITreeCluster<Double> nc = new STITreeCluster<Double>(cl);
					nc.setData(weight);
					clusters.add(nc);
				}
				else{
					STITreeCluster<Double> nc = clusters.get(pos);
					nc.setData(nc.getData()+weight);
				}
			}
		}


		//order the clusters by its frequency
		for(int i=1;i<clusters.size();i++){
			STITreeCluster<Double> cl1 = clusters.get(i);
			for(int j=0;j<i;j++){
				STITreeCluster<Double> cl2 = clusters.get(j);
				if(cl2.getData() < cl1.getData()){
					clusters.remove(cl1);
					clusters.add(j,cl1);
					break;
				}
			}
		}
	/*
		for(STITreeCluster<Integer> c: clusters){
			System.out.println(c);
		}
	 */
		//get the clusters for result trees
		List<STITreeCluster> minClusters = new ArrayList<STITreeCluster>();
        if(percentage==0){
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
        }
        else{
            for(STITreeCluster<Double> cl: clusters){
                if(cl.getData()/size > percentage/100.0){
                    minClusters.add(cl);
                }
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

    public Tree inferSpeciesTreeUnrooted(List<MutableTuple<Tree,Double>> trees, Map<String,String> taxonMap, int percentage){
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

		List<STITreeCluster<Double>> clusters = new LinkedList<STITreeCluster<Double>>();

        double size = 0;
        for(MutableTuple<Tree,Double> tuple: trees){
            Tree tr = tuple.Item1;
            double weight = tuple.Item2;
            size += weight;
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
					STITreeCluster<Double> cd = new STITreeCluster<Double>(stcl);
					cd.setData(frequency);
					clusters.add(cd);
				}
				else{
					STITreeCluster<Double> cd = clusters.get(pos);
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
			STITreeCluster<Double> cd = clusters.get(i);
			changed[i] = 1;
			if(cd.getClusterSize() == stTaxa.length-1){
				cd.setData(cd.getData()+num_trees);
			}
			else{

				int pos = clusters.indexOf(cd.complementaryCluster());
				if(pos != -1){
					STITreeCluster<Double> cc = clusters.get(pos);
					cd.setData(cd.getData()+cc.getData());
					cc.setData(cd.getData());
					changed[pos] = 1;
				}
			}
		}
*/

		//order the clusters by its frequency
		for(int i=1;i<clusters.size();i++){
			STITreeCluster<Double> cl1 = clusters.get(i);
			for(int j=0;j<i;j++){
				STITreeCluster<Double> cl2 = clusters.get(j);
				if(cl2.getData() < cl1.getData()){
					clusters.remove(cl1);
					clusters.add(j,cl1);
					break;
				}
			}
		}
/*
		for(STITreeCluster<Double> c: clusters){
			System.out.println(c);
		}
*/
		//get the clusters for result trees
		List<STITreeCluster> minClusters = new ArrayList<STITreeCluster>();

        if(percentage == 0){
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

        }
        else{
            for(STITreeCluster<Double> cl: clusters){
                if(cl.getData()/size> percentage/100.0){
                    minClusters.add(cl);
                }
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
