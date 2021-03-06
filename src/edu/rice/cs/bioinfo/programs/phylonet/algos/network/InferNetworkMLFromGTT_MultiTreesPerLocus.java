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
 * Created with Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This class is a subclass of InferNetworkMLFromGTT. It handles the cases each locus has multiple gene trees.
 */
public class InferNetworkMLFromGTT_MultiTreesPerLocus extends InferNetworkMLFromGTT {

    /**
     * Constructor of the class
     */
    public InferNetworkMLFromGTT_MultiTreesPerLocus(){
        _likelihoodCalculator = new NetworkLikelihoodFromGTT_MultiTreesPerLocus();
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
        int treeID = 0;
        Map<String, MutableTuple<Integer,Double>> tree2Info = new HashMap<String, MutableTuple<Integer,Double>>();
        for(Object o: originalGTs) {
            List<MutableTuple<Tree,Double>> treesForOneLocus = (List<MutableTuple<Tree,Double>>)o;
            Map<String, Integer> tree2infoIndex = new HashMap<String, Integer>();
            List<MutableTuple<Integer,Double>> infoList = new ArrayList<MutableTuple<Integer, Double>>();
            for (MutableTuple<Tree, Double> gtTuple : treesForOneLocus) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Integer,Double> existingInfo = tree2Info.get(exp);
                if (existingInfo == null) {
                    existingInfo = new MutableTuple<Integer,Double>(treeID, gtTuple.Item2);
                    gtsForStartingNetwork.add(gtTuple);
                    gtsForInferNetwork.add(gtTuple.Item1);
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
            MutableTuple<Integer,Double> tuple = (MutableTuple<Integer,Double>)gtsForStartingNetwork.get(info.Item1);
            tuple.Item2 = info.Item2;
        }
    }




}
