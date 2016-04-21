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

import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.HillClimberBase;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.SimpleHillClimbing;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing.SimulatedAnnealingSalterPearL;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NonUltrametricNetworkRandomParameterNeighbourGenerator;

import java.util.*;

/**
 * Created with Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is a subclass of InferNetworkMLFromGT. It's uses both topologies and branch lengths of gene trees as input data to infer species networks.
 * Note that gene trees need to be ultrametric and the branch lengths are in coalescent unit
 *
 * See "Maximum Likelihood Inference of Reticulate Evolutionary Histories", Proceedings of the National Academy of Sciences, 2014
 */
public abstract class InferNetworkMLFromGTTBL extends InferNetworkMLFromGT {

    /**
     * This function is to find the set of branches whose lengths cannot be estimated so that they can be ignored during the inference
     * In this case, none of the branches whose lengths can be ignored during the inference
     *
     * @param speciesNetwork        the species network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param singleAlleleSpecies   species that have only one allele sampled from it
     */
    protected void findSingleAlleleSpeciesSet(Network speciesNetwork, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies){}


    /**
     * This function is to obtain the random parameter generator for a species network
     *
     * @param dataForNetworkInference   (distinct) data used for network inference
     * @param allele2species            mapping from allele to the species it is sampled from
     * @param singleAlleleSpecies       species that have only one allele sampled from each
     */
    protected NetworkRandomParameterNeighbourGenerator getNetworkRandomParameterNeighbourGenerator(List dataForNetworkInference, Map<String,String> allele2species, Set<String> singleAlleleSpecies){
        return new NonUltrametricNetworkRandomParameterNeighbourGenerator();
    }


    /**
     * This function is to obtain the searching strategy which is simple hill climbing in this case where branch lengths are optimized for every network being tried
     * For using both branch lengths and topologies of gene trees, sampling branch lengths, which is used in simulated annealing, is not implemented
     */
    protected HillClimberBase getSearchingStrategy(Comparator<Double> comparator, NetworkRandomNeighbourGenerator allNeighboursStrategy){
        return new SimpleHillClimbing(comparator, allNeighboursStrategy);
    }
}
