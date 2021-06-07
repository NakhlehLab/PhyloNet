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

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_Rooted;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created with Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is a subclass of InferNetworkML. It's a version that uses gene trees as input data to infer species networks.
 *
 * See "Maximum Likelihood Inference of Reticulate Evolutionary Histories", Proceedings of the National Academy of Sciences, 2014
 */
public abstract class InferNetworkMLFromGT extends InferNetworkML {

    /**
     * This function is to obtain starting network for the search
     * MDC on trees is used by default
     *
     * @param gts                 gene trees for inferring the starting network
     * @param species2alleles     mapping from species to alleles sampled from it
     * @param hybridSpecies       species under reticulation nodes in the species network
     * @param startingNetwork     starting network if specified by the users
     *
     * @return  starting network
     */
    protected String getStartNetwork(List gts, Map<String,List<String>> species2alleles, Set<String> hybridSpecies, Network<Object> startingNetwork){
        Set<String> existingHybrids = new HashSet<>();
        if(startingNetwork == null){
            Map<String,String> allele2species = null;
            if(species2alleles!=null){
                allele2species = new HashMap<String, String>();
                for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                    String species = entry.getKey();
                    for(String allele: entry.getValue()){
                        allele2species.put(allele,species);
                    }
                }
            }
            MDCInference_Rooted mdc = new MDCInference_Rooted();
            Solution sol;
            if(allele2species==null){
                sol = (Solution)(mdc.inferSpeciesTree(gts, false, 1, false, true, -1).get(0));
            }
            else{
                sol = (Solution)(mdc.inferSpeciesTree(gts, allele2species, false, 1, false, true, -1).get(0));
            }
            Tree startingTree= Trees.generateRandomBinaryResolution(sol._st);
            startingNetwork = edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(startingTree.toString());
        }else{
            checkNetworkWithHybrids(startingNetwork, existingHybrids);
        }

        for(String hybrid: hybridSpecies){
            if(!existingHybrids.contains(hybrid))
                createHybrid(startingNetwork, hybrid);
        }

        for(NetNode<Object> node: startingNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                if(node.getParentDistance(parent) == NetNode.NO_DISTANCE)
                    node.setParentDistance(parent,1.0);
                if(node.isNetworkNode()){
                    if(node.getParentProbability(parent) == NetNode.NO_PROBABILITY)
                        node.setParentProbability(parent, 0.5);
                }
            }
        }

        return startingNetwork.toString();
    }






}
