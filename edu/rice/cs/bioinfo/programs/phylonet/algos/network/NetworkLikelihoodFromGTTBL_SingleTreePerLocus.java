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
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Created by Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is a subclass of NetworkLikelihoodFromGTTBL. It handles the cases each locus has one gene tree.
 */
public class NetworkLikelihoodFromGTTBL_SingleTreePerLocus extends NetworkLikelihoodFromGTTBL {


    /**
     * This function is to summarize the input gene trees
     *
     * @param originalGTs               original input gene trees
     * @param allele2species            mapping from allele to species which it is sampled from
     * @param gtsForInferNetwork        summarized gene trees
     * @param treeCorrespondences       relationships between the original gene trees and the gene trees in summarizedGTs
     */
    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List gtsForInferNetwork, List treeCorrespondences){
        for(Object list: originalGTs) {
            for (MutableTuple<Tree, Double> gtTuple : (List<MutableTuple<Tree, Double>>)list) {
                gtsForInferNetwork.add(gtTuple.Item1);
                treeCorrespondences.add(gtTuple.Item2);
            }
        }
    }



    /**
     * This function is to calculate the final log likelihood using the correspondences between the summarized gene trees and the original gene trees
     *
     * @param probList            the probabilities of each summarized gene tree respectively
     * @param gtCorrespondences   the correspondences between the summarized gene trees and the original gene trees
     */
    protected double calculateFinalLikelihood(double[] probList, List gtCorrespondences){
        double totalProb = 0;
        Iterator weight = gtCorrespondences.iterator();
        for(double prob: probList){
            totalProb += Math.log(prob) * (Double)weight.next();
        }
        return totalProb;
    }


}
