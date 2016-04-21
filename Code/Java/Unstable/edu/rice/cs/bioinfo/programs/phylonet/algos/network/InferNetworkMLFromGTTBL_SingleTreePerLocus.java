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
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created with Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is a subclass of InferNetworkMLFromGTTBL. It handles the cases each locus has single gene tree.
 */
public class InferNetworkMLFromGTTBL_SingleTreePerLocus extends InferNetworkMLFromGTTBL {

    /**
     * Constructor of the class
     */
    public InferNetworkMLFromGTTBL_SingleTreePerLocus(){
        _likelihoodCalculator = new NetworkLikelihoodFromGTTBL_SingleTreePerLocus();
    }


    /**
     * This function is to summarize the input gene trees
     * Since branch lengths of gene trees are used, every gene tree is distinct
     *
     * @param originalGTs               original input data
     * @param allele2species            mapping from allele to species which it is sampled from
     * @param gtsForStartingNetwork    data for inferring the starting network
     * @param gtsForInferNetwork       (distinct) data used during the search
     * @param treeCorrespondences       relationships between the original data and the data in dataForInferNetwork
     */
    protected void summarizeData(List originalGTs, Map<String, String> allele2species, List gtsForStartingNetwork, List gtsForInferNetwork, List treeCorrespondences) {
        for (Object list : originalGTs) {
            for (MutableTuple<Tree, Double> gtTuple : (List<MutableTuple<Tree, Double>>) list) {
                gtsForInferNetwork.add(gtTuple.Item1);
                treeCorrespondences.add(gtTuple.Item2);
            }
            gtsForStartingNetwork.addAll((List<MutableTuple<Tree, Double>>) list);
        }

    }


}
