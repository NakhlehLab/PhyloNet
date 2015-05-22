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

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkLikelihoodFromGTTBL_MultiTreesPerLocus extends NetworkLikelihoodFromGTTBL {


    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List gtsForInferNetwork, List treeCorrespondences){
        int treeID = 0;
        for(Object o: originalGTs) {
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            List<MutableTuple<Integer,Double>> infoList = new ArrayList<MutableTuple<Integer, Double>>();
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                gtsForInferNetwork.add(gtTuple.Item1);
                infoList.add(new MutableTuple<Integer, Double>(treeID, gtTuple.Item2));
                treeID++;
            }
            treeCorrespondences.add(infoList);
        }

    }


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
