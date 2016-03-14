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
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 12/2/11
 * Time: 1:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class MDCURInference_ILP {

      private String _outputPath = "";
	private static boolean _print = true;


	public Solution inferSpeciesTree(String gurobi_path,List<Tree> trees){
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

		List<Tree> allTrees = new ArrayList<Tree>();
		int[] treeNumVector = new int[trees.size()];
		for(int i=0;i<treeNumVector.length;i++){
			treeNumVector[i]=0;
		}
	    List<STITreeCluster<int[]>> allClusters = new ArrayList<STITreeCluster<int[]>>();
	    index=0;
		for(Tree tr : trees){
				for(TNode n:tr.postTraverse()){
					((STINode)n).setParentDistance(TMutableNode.NO_DISTANCE);
				}
				int tree_num = 0;
				for(Tree t: tr.getAllRootingTrees()){
					tree_num++;
					allTrees.add(t);
					for(STITreeCluster c:t.getClusters(taxa, true)){
						if(!allClusters.contains(c)){
							allClusters.add(new STITreeCluster<int[]>(c));
						}
					}
				}
				treeNumVector[index++] = tree_num;
			}


		int max_coal = 0;
		int totalTreeNum = allTrees.size();

		for(STITreeCluster<int[]> cl : allClusters){
			index = 0;
			int[] coalVector = new int[totalTreeNum];
			for(Tree tr: allTrees){
				int coal_num = DeepCoalescencesCounter.getClusterCoalNum_rooted(tr, cl);
				coalVector[index++] = coal_num;
				if(coal_num > max_coal){
					max_coal = coal_num;
				}
			}
			cl.setData(coalVector);
		}

		max_coal++;

		for(STITreeCluster<int[]> cl : allClusters){
			int[] coalVector = cl.getData();
			for(int i=0;i<coalVector.length;i++){
				coalVector[i] = max_coal - coalVector[i];
			}
			cl.setData(coalVector);
		}

		Solution sol = null;
		try{
			sol = solveQP(gurobi_path,treeNumVector,allClusters,allTrees);

		}catch(IOException e){
			System.err.println(e.getMessage());
			e.getStackTrace();
		}
		return sol;
	}

	public Solution inferSpeciesTree(String gurobi_path,List<Tree> trees, Map<String,String> taxonMap){
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

		List<Tree> allTrees = new ArrayList<Tree>();
		int[] treeNumVector = new int[trees.size()];
		for(int i=0;i<treeNumVector.length;i++){
			treeNumVector[i]=0;
		}
	    List<STITreeCluster<int[]>> allClusters = new ArrayList<STITreeCluster<int[]>>();

	    int index = 0;
		for(Tree tr : trees){
			((STITree)tr).setRooted(false);
			if(tr.isRooted()){
				allTrees.add(tr);
				treeNumVector[index++] = 1;
				for(STITreeCluster gtCluster:tr.getClusters(gtTaxa, true)){
					STITreeCluster<int[]> stCluster = new STITreeCluster<int[]>(stTaxa);
					for (String s : gtCluster.getClusterLeaves()) {
						stCluster.addLeaf(taxonMap.get(s));
					}
					if(stCluster.getClusterSize()>1 && stCluster.getClusterSize()<stTaxa.length){
						if(!allClusters.contains(stCluster)){
							allClusters.add(stCluster);
						}
					}
				}
			}
			else{
				for(TNode n:tr.postTraverse()){
					((STINode)n).setParentDistance(TMutableNode.NO_DISTANCE);
				}
				int tree_num = 0;
				for(Tree t: tr.getAllRootingTrees()){
					tree_num++;
					allTrees.add(t);
					for(STITreeCluster gtCluster:tr.getClusters(gtTaxa, true)){
						STITreeCluster<int[]> stCluster = new STITreeCluster<int[]>(stTaxa);
						for (String s : gtCluster.getClusterLeaves()) {
							stCluster.addLeaf(taxonMap.get(s));
						}
						if(stCluster.getClusterSize()>1 && stCluster.getClusterSize()<stTaxa.length){
							if(!allClusters.contains(stCluster)){
								allClusters.add(stCluster);
							}
						}
					}
				}
				treeNumVector[index++] = tree_num;
			}
		}

		int max_coal = 0;
		int totalTreeNum = allTrees.size();
		for(STITreeCluster<int[]> cl : allClusters){
			index = 0;
			int[] coalVector = new int[totalTreeNum];
			for(Tree tr: allTrees){
				int coal_num = DeepCoalescencesCounter.getClusterCoalNum_rooted(tr,cl,taxonMap);
				coalVector[index++] = coal_num;
				if(coal_num > max_coal){
					max_coal = coal_num;
				}
				cl.setData(coalVector);
			}
		}

		max_coal++;

		for(STITreeCluster<int[]> cl : allClusters){
			int[] coalVector = cl.getData();
			for(int i=0;i<coalVector.length;i++){
				coalVector[i] = max_coal - coalVector[i];
			}
			cl.setData(coalVector);
		}

		Solution sol = null;
		try{
			sol = solveQP(gurobi_path,treeNumVector,allClusters,allTrees);

		}catch(IOException e){
			System.err.println(e.getMessage());
			e.getStackTrace();
		}
		return sol;
	}


	public void setPath(String path){
		_outputPath = path;
	}

	private Solution solveQP (String gurobi_path,int[] treeNumVector, List<STITreeCluster<int[]>> allClusters, List<Tree> allTrees) throws IOException{

		File gurobiInput =  File.createTempFile("gurobiInput", ".lp");
        gurobiInput.deleteOnExit();
		File gurobiOutput = File.createTempFile("gurobiOutput", "txt");
        gurobiOutput.deleteOnExit();
		File gurobiResult = File.createTempFile("gurobiResult", "txt");
        gurobiResult.deleteOnExit();
		File gurobiError =  File.createTempFile("gurobiError", "txt");
        gurobiError.deleteOnExit();

        /*if(gurobiInput.exists()){
			gurobiInput.delete();
		}
		if(gurobiOutput.exists()){
			gurobiOutput.delete();
		}
		if(gurobiResult.exists()){
			gurobiResult.delete();
		}
		if(gurobiError.exists()){
			gurobiError.delete();
		}
		gurobiInput.createNewFile();
		gurobiOutput.createNewFile();
		gurobiResult.createNewFile(); */

		try {

			FileWriter fw = new FileWriter(gurobiInput);
			try {
				//System.out.println(generateGurobiInput(treeNumVector, allClusters).toString());
				fw.write(generateGurobiInput(treeNumVector, allClusters).toString());
			}
			finally {
				fw.close();
			}

			//String[] cmdArray = {"/bin/sh","-c",gurobi_path};
			//String[] cmdArray = {gurobi_path};
			//String cmdArray = gurobi_path;
			Process gurobi = Runtime.getRuntime().exec(gurobi_path);

			BufferedWriter gurobiCmdStream = new BufferedWriter(new OutputStreamWriter(gurobi.getOutputStream()));
			try {
				gurobiCmdStream.write("import gurobipy\n");
				gurobiCmdStream.flush();
				gurobiCmdStream.write("m = gurobipy.read(\""+gurobiInput.getPath()+"\")\n");
				//gurobiCmdStream.write("help()\n");
				gurobiCmdStream.flush();
				gurobiCmdStream.write("m.optimize()\n");
				gurobiCmdStream.flush();
				gurobiCmdStream.write("m.printAttr(\'X\')\n");
			} finally {
				gurobiCmdStream.close();
			}

			BufferedReader gurobiOutStream = new BufferedReader(new InputStreamReader(gurobi.getInputStream()));
			try {
				fw = new FileWriter(gurobiOutput);
				try {
					String line;
					while ((line = gurobiOutStream.readLine()) != null) {
						//System.out.println("output:"+line);
						fw.write(line + "\n");
					}
				} finally {
					fw.close();
				}
			} finally {
				gurobiOutStream.close();
			}


			BufferedReader cplexErrorStream = new BufferedReader(new InputStreamReader(gurobi.getErrorStream()));
			try {
				FileWriter fw2 = new FileWriter(gurobiError);
				try{
					String line;
					while ((line = cplexErrorStream.readLine()) != null) {
						fw2.write(line + "\n");
					}
				}finally{
					fw2.close();
				}
			} finally {
				cplexErrorStream.close();
			}

			try {
				if (gurobi.waitFor() != 0) {
					throw new RuntimeException("Abnormal termination of GUROBI.");
				}
			} catch (InterruptedException e) {
				throw new RuntimeException("GUROBI was interrupted during execution.", e);
			}

			BufferedReader br = new BufferedReader(new FileReader(gurobiOutput));
			String line;
			boolean start = false;
			List<STITreeCluster> minClusters = new ArrayList<STITreeCluster>();
			List<Integer> geneTrees = new ArrayList<Integer>();
			while ((line = br.readLine()) != null) {
				if(start){
					line = line.trim();
					int index = line.indexOf(" ");
					String variable = line.substring(0,index);
					int id = Integer.parseInt(variable.substring(1));
					if(variable.startsWith("x")){
						minClusters.add(allClusters.get(id));
					}
					else if(variable.startsWith("y")){
						geneTrees.add(id);
					}
				}
				if(line.startsWith("----")){
					start = true;
				}
			}
			Solution sol = new Solution();
			sol._st = Trees.buildTreeFromClusters(minClusters);
			List<MutableTuple<Tree,Double>> selectedgts = new ArrayList<MutableTuple<Tree,Double>>();
			for(int i:geneTrees){
				selectedgts.add(new MutableTuple<Tree, Double>(allTrees.get(i),1.0));
			}
			sol._totalCoals = DeepCoalescencesCounter.countExtraCoal(selectedgts, sol._st, true, 1);
			return sol;

		} finally {

			// Clean temporary files, if we got far enough to create them.

		/*	if (gurobiInput.exists() && !gurobiInput.delete()) {
				System.err.println("Cannot delete file " + gurobiInput + ". You may need to remove this file manually.");
			}
			if (gurobiOutput.exists() && !gurobiOutput.delete()) {
				System.err.println("Cannot delete file " + gurobiOutput + ". You may need to remove this file manually.");
			}    */

		}
	}

	/*
	public int countExtraCoal(List<Tree> gts,Tree st){
		int sum = 0;

		for(Tree gt : gts){
			sum += countExtraCoal(gt,st);
		}

		return sum;
	}

	public int countExtraCoal(Tree gt,Tree st){
		int min_coal = -1;

		List<Tree> gts = new ArrayList<Tree>();
		for(Tree tr: gt.getAllRootingTrees()){
			gts.clear();
			gts.add(tr);
			int coal = DeepCoalescencesCounter.countExtraCoal(gts, st, true);
			if(coal < min_coal || min_coal ==-1){
				min_coal = coal;
			}
		}

		return min_coal;
	}

	public int countExtraCoal(Tree gt,Tree st,Map<String,String> taxonMap){
		int min_coal = -1;

		List<Tree> gts = new ArrayList<Tree>();
		for(Tree tr: gt.getAllRootingTrees()){
			gts.clear();
			gts.add(tr);
			int coal = DeepCoalescencesCounter.countExtraCoal(gts, st,taxonMap, true);
			if(coal < min_coal || min_coal ==-1){
				min_coal = coal;
			}
		}

		return min_coal;
	}

	public int countExtraCoal(List<Tree> gts,Tree st,Map<String,String> taxonMap){
		int sum = 0;

		for(Tree gt : gts){
			sum += countExtraCoal(gt,st,taxonMap);
		}

		return sum;
	}
	*/

	private StringBuffer generateGurobiInput(int[] treeNumVector, List<STITreeCluster<int[]>> allClusters) {
		StringBuffer objective = new StringBuffer();
		StringBuffer constraint = new StringBuffer();
		StringBuffer variable = new StringBuffer();
		objective.append("Maximize\n ");
		constraint.append("Subject To\n");
		variable.append("Binaries\n");

		int count = 0;		// For line breaking.
		int i=0;
		int index = 0;
		int totalTree = 0;
		for (STITreeCluster<int[]> cl : allClusters) {
			int[] cv = cl.getData();
			for(int j=0;j<cv.length;j++){
				if (count > 0) {
					objective.append(" + ");
				}
				totalTree = cv.length;
				//objective.append(cv[j] + "z" + i+j);
				//constraint.append(" x"+i + "+y"+j + "-2z"+i+j + ">=0\n");
				//variable.append(" x"+i + " y"+j + " z"+i+j);
				objective.append(cv[j] + " z" + index);
				constraint.append(" x"+i + " + y"+j + " - 2 z"+index + " >= 0\n");
				//constraint.append(" x"+i + " + y"+j + " - z"+index + " <= 1\n");
				variable.append(" z"+index);
				index++;
				count++;
				if (count >= 50) {
					objective.append("\n");
					//variable.append("\n");
					count = 1;
				}
			}
			i++;
		}
		objective.append("\n");

		for(i=0;i<allClusters.size();i++){
			variable.append(" x"+i);
		}
		for(i=0;i<totalTree;i++){
			variable.append(" y"+i);
		}
		variable.append("\n");

		//generate more constraints
		//List<CompatibleClusterPair> compPairs = new LinkedList<CompatibleClusterPair>();
		for(i=0;i<allClusters.size();i++){
			STITreeCluster cl1 = allClusters.get(i);
			for(int j=i+1;j<allClusters.size();j++){
				STITreeCluster cl2 = allClusters.get(j);
				if(!cl1.isCompatible(cl2)){
					//compPairs.add(new CompatibleClusterPair(i,j));
					constraint.append(" x"+i + " + x"+j + " <= 1\n");
				}
			}
		}

		index = 0;
		for(i=0;i<treeNumVector.length;i++){
			int size = treeNumVector[i];
			constraint.append(" ");
			for(int j=0;j<size;j++){
				if(j>0){
					constraint.append(" + ");
				}
				constraint.append("y"+(index++));
			}
			constraint.append(" = 1\n");
		}

		StringBuffer formulation = new StringBuffer();
		formulation.append(objective);
		formulation.append(constraint);
		formulation.append(variable);
		formulation.append("End\n");

		return formulation;
	}


	class CompatibleClusterPair{
		int _cluster1;
		int _cluster2;

		public CompatibleClusterPair(int c1,int c2){
			_cluster1 = c1;
			_cluster2 = c2;
		}

		public boolean equals(Object o){
			if(!(o instanceof CompatibleClusterPair)){
				return false;
			}
			CompatibleClusterPair ccp = (CompatibleClusterPair)o;
			if ((ccp._cluster1 == _cluster1 && ccp._cluster2 == _cluster2)||
					(ccp._cluster2 == _cluster1 && ccp._cluster1 == _cluster2)) {
				return true;
			}
			else {
				return false;

			}
		}

		public int hashCode(){
			return (String.valueOf(_cluster1)).hashCode()+(String.valueOf(_cluster2)).hashCode();
		}
	}
}
