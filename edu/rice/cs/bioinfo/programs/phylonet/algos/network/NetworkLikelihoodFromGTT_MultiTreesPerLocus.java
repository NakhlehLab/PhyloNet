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
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is a subclass of NetworkLikelihoodFromGTT. It handles the cases each locus has multiple gene trees.
 */
public class NetworkLikelihoodFromGTT_MultiTreesPerLocus extends NetworkLikelihoodFromGTT {

    /**
     * This function is to summarize the input gene trees
     *
     * @param originalGTs               original input gene trees
     * @param allele2species            mapping from allele to species which it is sampled from
     * @param summarizedGTs             summarized gene trees
     * @param treeCorrespondences       relationships between the original gene trees and the gene trees in summarizedGTs
     */
    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List summarizedGTs, List treeCorrespondences){
        int treeID = 0;
        Map<String, Integer> tree2Info = new HashMap<String, Integer>();
        for(Object o: originalGTs) {
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            Map<String, Integer> tree2infoIndex = new HashMap<String, Integer>();
            List<MutableTuple<Integer,Double>> infoList = new ArrayList<MutableTuple<Integer, Double>>();
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                Integer existingID = tree2Info.get(exp);
                if (existingID == null) {
                    existingID = treeID;
                    summarizedGTs.add(gtTuple.Item1);
                    tree2Info.put(exp, existingID);
                    tree2infoIndex.put(exp, infoList.size());
                    infoList.add(new MutableTuple(treeID, gtTuple.Item2));
                    treeID++;
                } else {
                    Integer infoID = tree2infoIndex.get(exp);
                    if(infoID == null){
                        tree2infoIndex.put(exp, infoList.size());
                        infoList.add(new MutableTuple(existingID, gtTuple.Item2));
                    }
                    else{
                        infoList.get(infoID).Item2 += gtTuple.Item2;
                    }
                }



            }
            treeCorrespondences.add(infoList);
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
        for(Object o: gtCorrespondences){
            List<MutableTuple<Integer,Double>> locusInfo = (List<MutableTuple<Integer,Double>>)o;
            double totalProbForLocus = 0;
            double totalWeight = 0;
            for(MutableTuple<Integer,Double> one: locusInfo){
                totalProbForLocus += probList[one.Item1] * one.Item2;
                totalWeight += one.Item2;
            }
            totalProbForLocus = totalProbForLocus/totalWeight;
            totalProb += Math.log(totalProbForLocus);
        }
        return totalProb;
    }


}
