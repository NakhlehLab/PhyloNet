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
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.*;

/**
 * Created with Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is a subclass of InferNetworkMLFromGTT. It handles the cases each locus has single gene tree.
 */
public class InferNetworkMLFromGTT_SingleTreePerLocus extends InferNetworkMLFromGTT {

    /**
     * Constructor of the class
     */
    public InferNetworkMLFromGTT_SingleTreePerLocus(){
        _likelihoodCalculator = new NetworkLikelihoodFromGTT_SingleTreePerLocus();
    }


    /**
     * This function is to summarize the input gene trees by finding the distinct gene tree topologies
     *
     * @param originalGTs               original input data
     * @param allele2species            mapping from allele to species which it is sampled from
     * @param gtsForStartingNetwork    data for inferring the starting network
     * @param gtsForInferNetwork       (distinct) data used during the search
     * @param treeCorrespondences       relationships between the original data and the data in dataForInferNetwork
     */
    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List gtsForStartingNetwork, List gtsForInferNetwork, List treeCorrespondences){
        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        Map<String,Integer> exp2ID = new HashMap<String, Integer>();
        for(Object list: originalGTs) {
            for (MutableTuple<Tree, Double> gtTuple : (List<MutableTuple<Tree, Double>>)list) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
                if (existingTuple == null) {
                    existingTuple = new MutableTuple<>(gtTuple.Item1, gtTuple.Item2);
                    exp2tree.put(exp, existingTuple);
                    Set<Integer> binaryIDs = new HashSet<Integer>();
                    for (Tree btr : Trees.getAllBinaryResolution(gtTuple.Item1)) {
                        String btrExp = Trees.getLexicographicNewickString(btr, allele2species);
                        Integer index = exp2ID.get(btrExp);
                        if (index == null) {
                            index = gtsForInferNetwork.size();
                            gtsForInferNetwork.add(btr);
                            binaryIDs.add(index);
                            exp2ID.put(btrExp, index);
                        } else {
                            binaryIDs.add(index);
                        }
                    }
                    treeCorrespondences.add(new Tuple<MutableTuple<Tree, Double>, Set<Integer>>(existingTuple, binaryIDs));

                } else {
                    existingTuple.Item2 += gtTuple.Item2;
                }

            }
        }

        gtsForStartingNetwork.addAll(exp2tree.values());

    }


}
