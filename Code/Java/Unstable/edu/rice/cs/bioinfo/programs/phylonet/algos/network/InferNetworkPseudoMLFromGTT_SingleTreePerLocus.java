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

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created with Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is a subclass of InferNetworkPseudoMLFromGTT. It handles the cases each locus has single gene tree.
 */
public class InferNetworkPseudoMLFromGTT_SingleTreePerLocus extends InferNetworkPseudoMLFromGTT {

    /**
     * Constructor of the class
     */
    public InferNetworkPseudoMLFromGTT_SingleTreePerLocus(){
        _fullLikelihoodCalculator = new NetworkLikelihoodFromGTT_SingleTreePerLocus();
        _likelihoodCalculator = new NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus();
        if(_batchSize!=0){
            ((NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus)_likelihoodCalculator).setBatchSize(_batchSize);
        }
    }


    /**
     * This function is to summarize the input gene trees by finding the distinct gene tree topologies
     *
     * @param originalGTs               original input data
     * @param allele2species            mapping from allele to species which it is sampled from
     * @param gtsForStartingNetwork    data for inferring the starting network
     * @param allTriplets              all triplets
     * @param tripletFrequencies       triplet frequencies
     */
    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List gtsForStartingNetwork, List allTriplets, List tripletFrequencies){
        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        for(Object list: originalGTs) {
            for (MutableTuple<Tree, Double> gtTuple : (List<MutableTuple<Tree, Double>>)list) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
                if (existingTuple == null) {
                    existingTuple = gtTuple;
                    exp2tree.put(exp, existingTuple);

                } else {
                    existingTuple.Item2 += gtTuple.Item2;
                }

            }
        }
        gtsForStartingNetwork.addAll(exp2tree.values());
        NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus.computeTripleFrequenciesInGTs(gtsForStartingNetwork, allele2species, allTriplets, tripletFrequencies);
    }


}
