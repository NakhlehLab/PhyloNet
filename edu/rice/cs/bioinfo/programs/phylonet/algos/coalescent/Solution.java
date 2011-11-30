package edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

/**
 * This class saves the result of multiple versions of MDC
 *
 * @author Yun
 *
 */
public class Solution{
	public Tree _st;
	public int _totalCoals;
	int[] _clusterIDs;

	public Tree getTree(){
		return _st;
	}

	public int getCoalNum(){
		return _totalCoals;
	}
}


