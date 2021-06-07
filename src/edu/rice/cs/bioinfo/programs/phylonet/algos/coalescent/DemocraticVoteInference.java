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

import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by Yun Yu
 * Date: 11/14/11
 * Time: 5:07 PM
 *
 * This class is to use democratic vote to infer species tree from gene trees: take the gene tree with the highest frequency
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
		for (Tree tr : trees){
			boolean found=false;
			for (Map.Entry<Tree, Counter> entry : treeCounter.entrySet()){
				if(Trees.haveSameRootedTopology(entry.getKey(), tr)){
					entry.getValue().increase();
					found = true;
					break;
				}
			}
			if(!found){
				treeCounter.put(tr, new Counter());
			}
		}
		return treeCounter;
	}


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
