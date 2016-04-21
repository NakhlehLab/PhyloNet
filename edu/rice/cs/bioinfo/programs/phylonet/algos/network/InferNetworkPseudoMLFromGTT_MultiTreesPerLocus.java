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
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing.SimulatedAnnealingSalterPearL;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomTopologyNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created with Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is a subclass of InferNetworkPseudoMLFromGTT. It handles the cases each locus has multiple gene trees.
 */
public class InferNetworkPseudoMLFromGTT_MultiTreesPerLocus extends InferNetworkPseudoMLFromGTT {

    /**
     * Constructor of the class
     */
    public InferNetworkPseudoMLFromGTT_MultiTreesPerLocus(){
        _fullLikelihoodCalculator = new NetworkLikelihoodFromGTT_MultiTreesPerLocus();
        _likelihoodCalculator = new NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus();
        if(_batchSize!=0){
            ((NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus)_likelihoodCalculator).setBatchSize(_batchSize);
        }
    }


    /**
     * This function is to summarize the input gene trees by finding the distinct gene tree topologies
     *
     * @param originalGTs               original input data
     * @param allele2species            mapping from allele to species which it is sampled from
     * @param gtsForStartingNetwork     data for inferring the starting network
     * @param allTriplets               all triplets
     * @param tripletFrequencies        triplet frequencies
     */
    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List gtsForStartingNetwork, List allTriplets, List tripletFrequencies){
        int treeID = 0;
        Map<String, Integer> distinctTree2ID = new HashMap<>();
        List<List<MutableTuple<Integer,Double>>> treeCorrespondences = new ArrayList<>();
        int i=0;
        for(Object o: originalGTs) {
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            Map<String, Integer> tree2infoIndex = new HashMap<String, Integer>();
            List<MutableTuple<Integer,Double>> infoList = new ArrayList<MutableTuple<Integer, Double>>();
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                Integer existingTreeID = distinctTree2ID.get(exp);
                if (existingTreeID == null) {
                    existingTreeID = treeID;
                    gtsForStartingNetwork.add(new MutableTuple<Tree, Double>(gtTuple.Item1,gtTuple.Item2));
                    distinctTree2ID.put(exp, existingTreeID);
                    tree2infoIndex.put(exp, infoList.size());
                    infoList.add(new MutableTuple(treeID, gtTuple.Item2));
                    treeID++;
                } else {
                    ((MutableTuple<Tree,Double>)(gtsForStartingNetwork.get(existingTreeID))).Item2 += gtTuple.Item2;
                    Integer infoID = tree2infoIndex.get(exp);
                    if(infoID == null){
                        tree2infoIndex.put(exp, infoList.size());
                        infoList.add(new MutableTuple(existingTreeID, gtTuple.Item2));
                    }
                    else{
                        infoList.get(infoID).Item2 += gtTuple.Item2;
                    }
                }
            }
            treeCorrespondences.add(infoList);
        }
        NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus.computeTripleFrequenciesInGTs(gtsForStartingNetwork, allele2species, treeCorrespondences, allTriplets, tripletFrequencies);
    }


}
