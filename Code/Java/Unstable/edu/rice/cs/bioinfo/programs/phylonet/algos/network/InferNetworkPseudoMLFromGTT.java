/*
 * Copyright (c) 2013 Rice University.
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

package edu.rice.cs.bioinfo.programs.phylonet.algos.network;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing.SimulatedAnnealingSalterPearL;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomTopologyNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class InferNetworkPseudoMLFromGTT extends InferNetworkMLFromGTT {
    NetworkLikelihood _fullLikelihoodCalculator;
    protected int _batchSize;

    public void setBatchSize(int size) {
        _batchSize = size;
    }

    public void inferNetwork(List originalData, Map<String,List<String>> species2alleles, int maxReticulations, int numSol, boolean optimization, LinkedList<Tuple<Network,Double>> resultList){
        super.inferNetwork(originalData, species2alleles, maxReticulations, numSol, !optimization, resultList);
        if(optimization){
            Map<String,String> allele2species = null;
            if(species2alleles!=null){
                allele2species = new HashMap<String, String>();
                for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                    for(String allele: entry.getValue()){
                        allele2species.put(allele, entry.getKey());
                    }
                }
            }
            List summarizedGTs = new ArrayList();
            List gtCorrespondence = new ArrayList();
            _fullLikelihoodCalculator.summarizeData(originalData, allele2species, summarizedGTs, gtCorrespondence);
            Set<String> singleAlleleSpecies = new HashSet<>();
            findSingleAlleleSpeciesSet(resultList.get(0).Item1, singleAlleleSpecies);

            LinkedList<Tuple<Network,Double>> updatedResult = optimizeResultingNetworks(_fullLikelihoodCalculator, summarizedGTs, gtCorrespondence, species2alleles, singleAlleleSpecies, resultList);
            resultList.clear();
            resultList.addAll(updatedResult);
        }
    }

    private void findSingleAlleleSpeciesSet(Network network, Set<String> singleAlleleSpecies){
        for(Object node: network.getLeaves()){
            singleAlleleSpecies.add(((NetNode)node).getName());
        }
    }
}
