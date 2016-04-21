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
 * Created with Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is a subclass of InferNetworkMLFromGTT.
 * It infers species networks from a collection of gene trees under maximum pseudo-likelihood
 *
 * See "A Maximum Pseudo-likelihood Approach for Phylogenetic Networks", BMC Genomics, 2015.
 */
public abstract class InferNetworkPseudoMLFromGTT extends InferNetworkMLFromGTT {
    NetworkLikelihood _fullLikelihoodCalculator;
    int _batchSize;


    /**
     * This function is to set the batch size for computing the pseudo-likelihood in parallel
     */
    public void setBatchSize(int size) {
        _batchSize = size;
    }


    /**
     * This is the main function for inferring a species network from gene trees under maximum pseudo-likelihood
     *
     * @param gts                   the gene trees used for inferring the species network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param maxReticulations      the maximal number of reticulations in the inferred species network
     * @param numSol                number of solutions requested to return
     * @param postOptimization      whether optimizing branch lengths and inheritance probabilities of the inferred species networks under maximum FULL likelihood is needed
     *                              by default, the method optimizes branch lengths and inheritance probabilities of the inferred species networks under maximum pseudo-likelihood
     * @param resultList            resulting species networks along with their log likelihood
     */
    public void inferNetwork(List gts, Map<String,List<String>> species2alleles, int maxReticulations, int numSol, boolean postOptimization, LinkedList<Tuple<Network,Double>> resultList){
        super.inferNetwork(gts, species2alleles, maxReticulations, numSol, !postOptimization, resultList);
        if(postOptimization){
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
            _fullLikelihoodCalculator.summarizeData(gts, allele2species, summarizedGTs, gtCorrespondence);
            Set<String> singleAlleleSpecies = new HashSet<>();
            findSingleAlleleSpeciesSet(resultList.get(0).Item1, singleAlleleSpecies);

            LinkedList<Tuple<Network,Double>> updatedResult = optimizeResultingNetworks(_fullLikelihoodCalculator, summarizedGTs, gtCorrespondence, species2alleles, singleAlleleSpecies, resultList);
            resultList.clear();
            resultList.addAll(updatedResult);
        }
    }


    /**
     * This function is to find the set of branches whose lengths cannot be estimated so that they can be ignored during the inference
     * In this case, branches incident with nodes who has only one leaf node under it do not have any information about the lengths, so they can be ignored during the inference
     *
     * @param network               the species network
     * @param singleAlleleSpecies   species that have only one allele sampled from it
     */
    private void findSingleAlleleSpeciesSet(Network network, Set<String> singleAlleleSpecies){
        for(Object node: network.getLeaves()){
            singleAlleleSpecies.add(((NetNode)node).getName());
        }
    }
}
