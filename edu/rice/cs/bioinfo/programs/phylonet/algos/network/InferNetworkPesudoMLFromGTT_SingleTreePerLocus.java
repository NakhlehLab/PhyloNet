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
public class InferNetworkPesudoMLFromGTT_SingleTreePerLocus extends InferNetworkMLFromGTT {

    public void setNumRuns(int numRuns){
        _numRuns = numRuns;
    }


    protected void summarizeData(List originalGTs, Map<String,String> allele2species, List gtsForStartingNetwork, List tripleFrequencies, List treeCorrespondences){
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
        tripleFrequencies.addAll(computeTripleFrequenciesInGTs(gtsForStartingNetwork, allele2species));
    }


    private List<MutableTuple<String, double[]>> computeTripleFrequenciesInGTs(List<MutableTuple<Tree,Double>> gts, Map<String,String> allele2species){
        String[] speciesArray;
        if(allele2species==null){
            speciesArray = gts.get(0).Item1.getLeaves();
        }
        else{
            Set<String> taxonSet = new HashSet<>();
            taxonSet.addAll(allele2species.values());
            speciesArray = new String[taxonSet.size()];
            int index = 0;
            for(String taxon: taxonSet){
                speciesArray[index++] = taxon;
            }
        }
        Map<String, double[]> triple2counts = new HashMap<>();
        for(int i=0; i<speciesArray.length; i++){
            for(int j=i+1; j<speciesArray.length; j++){
                for(int k=j+1; k<speciesArray.length; k++){
                    triple2counts.put(speciesArray[i]+"&"+speciesArray[j]+"&"+speciesArray[k], new double[3]);
                }
            }
        }
        Map<String, Integer> species2ID = new HashMap<>();
        for(int i=0; i<speciesArray.length; i++){
            species2ID.put(speciesArray[i], i);
        }

        for(MutableTuple<Tree,Double> gt: gts){
            if(allele2species==null) {

                computeTripleFrequenciesFromSingleGT(gt, speciesArray, species2ID, triple2counts);
            }else{
                Map<String, Integer> allele2speciesID = new HashMap<>();
                for(Map.Entry<String,String> entry: allele2species.entrySet()){
                    allele2speciesID.put(entry.getKey(), species2ID.get(entry.getValue()));
                }
                computeTripleFrequenciesFromSingleGT(gt, allele2speciesID, speciesArray, triple2counts);
            }
        }
        List<MutableTuple<String, double[]>> tripleFrequencies = new ArrayList<>();
        for(Map.Entry<String, double[]> entry: triple2counts.entrySet()){
            tripleFrequencies.add(new MutableTuple<String, double[]>(entry.getKey(), entry.getValue()));
        }
        return tripleFrequencies;
    }



    private void computeTripleFrequenciesFromSingleGT(MutableTuple<Tree,Double> gt, String[] taxa, Map<String, Integer> taxon2ID, Map<String, double[]> triple2counts){
        int numTaxa = taxon2ID.size();
        int[][] pairwiseDepths = new int[numTaxa][numTaxa];
        Map<TNode, Integer> node2depth = new HashMap<>();
        Map<TNode, List<String>> node2leaves = new HashMap<>();
        for(TNode node: gt.Item1.postTraverse()){
            int depth = 0;
            List<String> leaves = new ArrayList<>();
            if(node.isLeaf()){
                //leaves.add(taxon2ID.get(node.getName()));
                leaves.add(node.getName());
            }
            else{
                List<List<String>> childLeavesList = new ArrayList<>();
                for(TNode child: node.getChildren()){
                    depth = Math.max(depth, node2depth.get(child));
                    List<String> childLeaves = new ArrayList<>();
                    childLeaves.addAll(node2leaves.get(child));
                    childLeavesList.add(childLeaves);
                    leaves.addAll(childLeaves);
                }
                depth++;
                for(int i=0; i<childLeavesList.size(); i++){
                    List<String> childLeaves1 = childLeavesList.get(i);
                    for(int j=i+1; j<childLeavesList.size(); j++){
                        List<String> childLeaves2 = childLeavesList.get(j);
                        for(String leaf1: childLeaves1){
                            int id1 = taxon2ID.get(leaf1);
                            for(String leaf2: childLeaves2){
                                int id2 = taxon2ID.get(leaf2);
                                pairwiseDepths[id1][id2] = depth;
                                pairwiseDepths[id2][id1] = depth;
                            }
                        }
                    }

                }

            }
            node2depth.put(node, depth);
            node2leaves.put(node, leaves);
        }

        for(int i=0; i<numTaxa; i++){
            for(int j=i+1; j<numTaxa; j++){
                int ij = pairwiseDepths[i][j];
                String pair = taxa[i]+"&"+taxa[j];
                for(int k=j+1; k<numTaxa; k++){
                    int minIndex = -1;
                    int ik = pairwiseDepths[i][k];
                    int jk = pairwiseDepths[j][k];
                    if(ij<ik && ij<jk){
                        minIndex = 0;
                    }else if(ik<ij && ik<jk)
                    {
                        minIndex = 1;
                    }
                    else if(jk<ij && jk<ik){
                        minIndex = 2;
                    }
                    if(minIndex!=-1){
                        triple2counts.get(pair+"&"+taxa[k])[minIndex]+=gt.Item2;
                    }
                    else{
                        double[] frequencies = triple2counts.get(pair+"&"+taxa[k]);
                        for(int m=0; m<3; m++){
                            frequencies[m] += gt.Item2/3;
                        }
                    }
                }
            }
        }
    }


    private void computeTripleFrequenciesFromSingleGT(MutableTuple<Tree,Double> gt, Map<String,Integer> allele2speciesID, String[] speciesArray, Map<String, double[]> triple2counts){
        Set<String> allAlleles = new HashSet<>();
        int[] alleleNums = new int[speciesArray.length];
        for(String allele: gt.Item1.getLeaves()){
            allAlleles.add(allele);
            alleleNums[allele2speciesID.get(allele)]++;
        }
        Map<TNode,Set<String>> node2leaves = new HashMap<>();
        for(TNode node: gt.Item1.postTraverse()){
            Set<String> leavesUnder = new HashSet<>();
            node2leaves.put(node, leavesUnder);
            if(node.isLeaf()){
                leavesUnder.add(node.getName());
            }
            else{
                List<Set<String>> childLeavesList = new ArrayList<>();
                for(TNode child: node.getChildren()){
                    Set<String> childLeaves = node2leaves.get(child);
                    leavesUnder.addAll(childLeaves);
                    childLeavesList.add(childLeaves);
                }

                allAlleles.removeAll(leavesUnder);

                for(int i=0; i<childLeavesList.size(); i++){
                    Set<String> childLeaves1 = childLeavesList.get(i);
                    for(int j=i+1; j<childLeavesList.size(); j++){
                        Set<String> childLeaves2 = childLeavesList.get(j);
                        for(String allele1: childLeaves1){
                            int species1 = allele2speciesID.get(allele1);
                            for(String allele2: childLeaves2){
                                int species2 = allele2speciesID.get(allele2);
                                if(species1!=species2){
                                    for(String allele3: allAlleles){
                                        int species3 = allele2speciesID.get(allele3);
                                        if(species1!=species3 && species2!=species3){
                                            addHighestFrequency(species1, species2, species3, speciesArray, alleleNums, triple2counts);
                                        }
                                    }
                                }
                            }
                        }
                        //non-binary node
                        for(int k=j+1; k<childLeavesList.size(); k++) {
                            Set<String> childLeaves3 = childLeavesList.get(k);
                            for(String allele1: childLeaves1) {
                                int species1 = allele2speciesID.get(allele1);
                                for (String allele2 : childLeaves2) {
                                    int species2 = allele2speciesID.get(allele2);
                                    for (String allele3 : childLeaves3) {
                                        int species3 = allele2speciesID.get(allele3);
                                        addEqualFrequency(species1, species2, species3, speciesArray, alleleNums, triple2counts);
                                    }
                                }
                            }
                        }
                    }

                }

                allAlleles.addAll(leavesUnder);
            }

        }

    }

    private void addEqualFrequency(int species1, int species2, int species3, String[] speciesArray, int[] alleleNums, Map<String, double[]> triple2counts) {
        double weight = 1.0/(alleleNums[species1]*alleleNums[species2]*alleleNums[species3]);
        int[] forSort = new int[3];
        forSort[0] = species1;
        forSort[1] = species2;
        forSort[2] = species3;
        Arrays.sort(forSort);
        String exp = speciesArray[forSort[0]]+"&"+speciesArray[forSort[1]]+"&"+speciesArray[forSort[2]];
        double[] frequency = triple2counts.get(exp);
        for(int i=0; i<3; i++){
            frequency[i] += weight/3;
        }
    }



    private void addHighestFrequency(int species1, int species2, int species3, String[] speciesArray, int[] alleleNums, Map<String, double[]> triple2counts) {
        double weight = 1.0/(alleleNums[species1]*alleleNums[species2]*alleleNums[species3]);
        if(species1<species3 && species2<species3){
            String exp = speciesArray[Math.min(species1,species2)]+"&"+speciesArray[Math.max(species1, species2)]+"&"+speciesArray[species3];
            triple2counts.get(exp)[0] += weight;
        }
        else if((species1<species3&&species3<species2) || (species2<species3&&species3<species1)){
            String exp = speciesArray[Math.min(species1,species2)]+"&"+speciesArray[species3]+"&"+speciesArray[Math.max(species1,species2)];
            triple2counts.get(exp)[1] += weight;
        }
        else if(species1>species3 && species2>species3){
            String exp = speciesArray[species3]+"&"+speciesArray[Math.min(species1,species2)]+"&"+speciesArray[Math.max(species1,species2)];
            triple2counts.get(exp)[2] += weight;
        }
    }




    protected double computeLikelihood(final Network<Object> speciesNetwork, final Map<String, List<String>> species2alleles, final List summarizedGTs, final List treeCorrespondence){
        //System.out.println(speciesNetwork.toString());
        NetworkPseudoLikelihoodFromGTT likelihoodComputer = new NetworkPseudoLikelihoodFromGTT();
        likelihoodComputer.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _numThreads);
        double prob = likelihoodComputer.computeLikelihood(speciesNetwork, species2alleles, summarizedGTs, treeCorrespondence, _optimizeBL);
        return prob;
    }



}
