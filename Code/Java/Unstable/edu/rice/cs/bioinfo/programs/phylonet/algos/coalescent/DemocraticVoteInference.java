package edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent;

import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/14/11
 * Time: 5:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class DemocraticVoteInference {

    private long _maxCount = 0;
	private static boolean _print = true;


	public long getFrequency(){
		return _maxCount;
	}

	public static void setPrinting(boolean print){
		_print = print;
	}

	public List<Tree> inferSpeciesTree(List<Tree> trees){
		Map<Tree,Counter> treeCounter = countGeneTrees(trees);
		List<Tree> inferredTrees = getMaxTrees(treeCounter);
	    _maxCount = Math.round((double)_maxCount*100/trees.size());
	    return inferredTrees;
	}

	private Map<Tree,Counter> countGeneTrees(List<Tree> trees){
		Map<Tree,Counter> treeCounter = new HashMap<Tree,Counter>();

		SymmetricDifference sd = new SymmetricDifference();

		for (Tree tr : trees){
			boolean found=false;
			for (Tree treeKey : treeCounter.keySet()){
				sd.computeDifference(tr, treeKey);
				if(sd.getFalseNegativeCount()==0 && sd.getFalsePositiveCount()==0){
					treeCounter.get(treeKey).increase();
					found=true;
					break;
				}
			}
			if(!found){
				treeCounter.put(tr, new Counter());
			}
		}
		return treeCounter;
		//System.out.println("mapsize:"+_treeCounter.size());
	}

	//private List<Tree> getMaxTrees(List<Tree> _inferredTrees){
		//List<Tree> _inferredTrees = new LinkedList<Tree>();
	private List<Tree> getMaxTrees(Map<Tree,Counter> treeCounter){
		int maxCount = 0;
		List<Tree> inferredTrees = new LinkedList<Tree>();

		for (Tree tr : treeCounter.keySet()){
			int count = treeCounter.get(tr).getCount();
			if(count > maxCount){
				inferredTrees.clear();
				inferredTrees.add(tr);
				maxCount = count;
			}
			else if(count == maxCount){
				inferredTrees.add(tr);
			}
		}
		_maxCount = maxCount;
		return inferredTrees;
	}

	class Counter{
		private int _num;

		public Counter(){
			_num = 1;
		}

		public void increase(){
			_num ++;
		}

		public int getCount(){
			return _num;
		}
	}
}
