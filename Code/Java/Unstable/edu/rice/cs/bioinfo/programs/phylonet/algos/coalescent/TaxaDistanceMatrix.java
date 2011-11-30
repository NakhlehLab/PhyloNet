package edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/14/11
 * Time: 1:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class TaxaDistanceMatrix
{
    private String[] taxa;
	private Map<STITreeCluster, Double> distanceMap;

	public TaxaDistanceMatrix(String[] taxa) {
		this.taxa = taxa;
		distanceMap = new HashMap<STITreeCluster, Double>();
	}

	public void put(STITreeCluster c, double d) {
		distanceMap.put(c, d);
	}

	public String[] getTaxa() {
		return taxa;
	}

	public Map<STITreeCluster, Double> getDistanceMap() {
		return distanceMap;
	}
}
