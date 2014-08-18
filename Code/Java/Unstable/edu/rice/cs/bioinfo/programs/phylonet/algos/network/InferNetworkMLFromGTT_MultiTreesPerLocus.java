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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferNetworkMLFromGTT_MultiTreesPerLocus extends InferNetworkMLFromGT {


    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List dataForStartingNetwork, List dataForInferNetwork, List treeCorrespondences){
        int treeID = 0;
        Map<String, MutableTuple<Integer,Double>> tree2Info = new HashMap<>();
        for(Object o: originalGTs) {
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            Map<String, Integer> tree2infoIndex = new HashMap<>();
            List<MutableTuple<Integer,Double>> infoList = new ArrayList<MutableTuple<Integer, Double>>();
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                //MutableTuple<Integer,Double> info = new MutableTuple<Integer, Double>(-1,gtTuple.Item2);
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Integer,Double> existingInfo = tree2Info.get(exp);
                if (existingInfo == null) {
                    existingInfo = new MutableTuple<>(treeID, gtTuple.Item2);
                    dataForStartingNetwork.add(gtTuple);
                    dataForInferNetwork.add(gtTuple.Item1);
                    tree2Info.put(exp, existingInfo);
                    tree2infoIndex.put(exp, infoList.size());
                    infoList.add(new MutableTuple(treeID, gtTuple.Item2));
                    treeID++;
                } else {
                    existingInfo.Item2 += gtTuple.Item2;
                    Integer infoID = tree2infoIndex.get(exp);
                    if(infoID == null){
                        tree2infoIndex.put(exp, infoList.size());
                        infoList.add(new MutableTuple(existingInfo.Item1, gtTuple.Item2));
                    }
                    else{
                        infoList.get(infoID).Item2 += gtTuple.Item2;
                    }
                }



            }
            treeCorrespondences.add(infoList);
        }

        for(MutableTuple<Integer,Double> info: tree2Info.values()){
            MutableTuple<Integer,Double> tuple = (MutableTuple<Integer,Double>)dataForStartingNetwork.get(info.Item1);
            tuple.Item2 = info.Item2;
        }

    }


    protected double findOptimalBranchLength(final Network<Object> speciesNetwork, final Map<String, List<String>> species2alleles, final List summarizedGTs, final List treeCorrespondence){
        NetworkLikelihoodFromGTT_MultiTreesPerLocus likelihoodComputer = new NetworkLikelihoodFromGTT_MultiTreesPerLocus();
        likelihoodComputer.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _numThreads);
        return likelihoodComputer.findOptimalBranchLength(speciesNetwork, species2alleles, summarizedGTs, treeCorrespondence);
    }



}