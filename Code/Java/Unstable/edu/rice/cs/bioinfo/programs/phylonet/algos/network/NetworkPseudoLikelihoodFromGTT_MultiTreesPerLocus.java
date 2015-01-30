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
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkPseudoLikelihoodFromGTT_MultiTreesPerLocus extends NetworkPseudoLikelihoodFromGTT {



    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List allTriplets, List tripletFrequencies){
        int treeID = 0;
        Map<String, Integer> distinctTree2ID = new HashMap<>();
        List<Tree> distinctGTs = new ArrayList<>();
        List<List<MutableTuple<Integer,Double>>> treeCorrespondences = new ArrayList<>();
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
                    distinctGTs.add(gtTuple.Item1);
                    distinctTree2ID.put(exp, existingTreeID);
                    tree2infoIndex.put(exp, infoList.size());
                    infoList.add(new MutableTuple(treeID, gtTuple.Item2));
                    treeID++;
                } else {
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
        computeTripleFrequenciesInGTs(distinctGTs, allele2species, treeCorrespondences, allTriplets, tripletFrequencies);
        System.out.println("Done");
    }


    private void computeTripleFrequenciesInGTs(List<Tree> gts, Map<String,String> allele2species, List<List<MutableTuple<Integer,Double>>> treeCorrespondences, List<String> allTriplets, List<List<Tuple<double[][],Double>>> tripletFrequencies){

        List<Map<String, double[]>> triple2countList = new ArrayList<>();
        Set<String> allAlleleSet = new HashSet<>();

        for(Tree gt: gts){
            Map<String, double[]> result;

            if(allele2species==null) {
                result = computeTripleFrequenciesFromSingleGT(gt);

            }else{
                result = computeTripleFrequenciesFromSingleGT(gt, allele2species);
            }
            allAlleleSet.addAll(result.keySet());
            triple2countList.add(result);
        }

        allTriplets.addAll(allAlleleSet);
        for(List<MutableTuple<Integer,Double>> oneLocus: treeCorrespondences){
            List<Tuple<double[][],Double>> freqList = new ArrayList<>();
            tripletFrequencies.add(freqList);
            for(MutableTuple<Integer,Double> gtInfo: oneLocus){
                double[][] gtTripleFreq = new double[allTriplets.size()][3];
                freqList.add(new Tuple<double[][], Double>(gtTripleFreq, gtInfo.Item2));
                Map<String, double[]> triple2count = triple2countList.get(gtInfo.Item1);
                int index = 0;
                for(String triplet: allTriplets){
                    double[] baseFreq = triple2count.get(triplet);
                    for(int i=0; i<3; i++){
                        gtTripleFreq[index][i] = baseFreq[i];
                    }
                    index++;
                }
            }
        }
    }




    protected double calculateFinalLikelihood(double[][] probs, List tripletFrequencies){

        double totalProb = 0;
        for(Object o: tripletFrequencies){
            List<Tuple<double[][],Double>> freqForOneLocus = (List<Tuple<double[][],Double>>)o;
            double probForOneLocus = 0;
            double totalWeight = 0;
            for(Tuple<double[][],Double> freqs: freqForOneLocus){
                double probForOneTree = freqs.Item2;
                totalWeight += freqs.Item2;
                for(int i=0; i<probs.length; i++){
                    for(int j=0; j<3; j++){
                        probForOneTree *= Math.pow(probs[i][j],freqs.Item1[i][j]);
                    }
                }
                probForOneLocus += probForOneTree / totalWeight;
            }
            totalProb += Math.log(probForOneLocus);

        }
        return totalProb;
    }

}
