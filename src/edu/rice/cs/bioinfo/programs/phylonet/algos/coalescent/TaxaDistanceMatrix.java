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
