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


import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.SimpleHillClimbing;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing.SimulatedAnnealingSalterPearL;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomTopologyNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NonUltrametricNetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferNetworkMLFromGTTWithCrossValidation extends InferNetworkMLFromGTT_SingleTreePerLocus {
    private Network[] _optimalNetworkswithReticulations;
    private double[] _optimalScoreswithReticulations;
    private int _numFolds;
    private boolean _printDetails = false;



    protected Func1<Network, Double> getScoreFunction(final List summarizedData, final Map<String, List<String>> species2alleles, final List dataCorrespondences, final Set<String> singleAlleleSpecies ){
        return new Func1<Network, Double>() {
            public Double execute(Network speciesNetwork) {

                double score = computeLikelihood(speciesNetwork, species2alleles, summarizedData, dataCorrespondences, singleAlleleSpecies);

                //System.out.println(speciesNetwork + ": " + score);
                // A score can be lower than the scores that I currently have in the
                // optimal score list, but it is still the best score
                // maybe for a lower reticulation network. Therefore, I just
                // compare every network found with those in the _optimalNetworkwithReticulautions
                // to update it with the best network of k reticulations found so far.

                // Save the best network with k reticulations into _optimalNetworkswithReticulations

                // Find the number of network nodes in speciesNetwork.

                int indexCount = speciesNetwork.getReticulationCount();

                if (score > _optimalScoreswithReticulations[indexCount]) {  // initially -inf, already assigned double value
                    _optimalNetworkswithReticulations[indexCount] = speciesNetwork.clone();
                    _optimalScoreswithReticulations[indexCount] = score;
                }
                System.gc();
                return score;
            }
        };
    }



    public void inferNetwork(List originalGTs, Map<String,List<String>> species2alleles, int maxReticulations, int numFolds, LinkedList<Tuple<Network,Double>> resultList){
        _numFolds = numFolds;
        _optimalNetworkswithReticulations = new Network[10];   // hold the best network so far with k reticulation nodes
        _optimalScoreswithReticulations = new double[10];
        Arrays.fill(_optimalScoreswithReticulations, Double.NEGATIVE_INFINITY);


        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new HashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                for(String allele: entry.getValue()){
                    allele2species.put(allele, entry.getKey());
                }
            }
        }

        List gtsForNetworkInference = new ArrayList();
        List gtsForStartingNetwork = new ArrayList();
        List gtCorrespondence = new ArrayList();
        summarizeData(originalGTs, allele2species, gtsForStartingNetwork, gtsForNetworkInference, gtCorrespondence);

        String startingNetwork = getStartNetwork(gtsForStartingNetwork, species2alleles, _fixedHybrid, _startNetwork);
        Network speciesNetwork = Networks.readNetwork(startingNetwork);
        gtsForStartingNetwork.clear();

        Set<String> singleAlleleSpecies = new HashSet<>();
        findSingleAlleleSpeciesSet(speciesNetwork, species2alleles, singleAlleleSpecies);

        NetworkRandomNeighbourGenerator allNeighboursStrategy = new NetworkRandomNeighbourGenerator(new NetworkRandomTopologyNeighbourGenerator(_topologyOperationWeight, maxReticulations, _moveDiameter, _reticulationDiameter, _fixedHybrid, _seed), _topologyVsParameterOperation[0], new NonUltrametricNetworkRandomParameterNeighbourGenerator(singleAlleleSpecies), _topologyVsParameterOperation[1], _seed);
        Comparator<Double> comparator = getDoubleScoreComparator();
        //SimpleHillClimbing searcher = new SimpleHillClimbing(comparator, allNeighboursStrategy);

        SimulatedAnnealingSalterPearL searcher = new SimulatedAnnealingSalterPearL(comparator, allNeighboursStrategy, _seed);
        searcher.setLogFile(_logFile);
        //System.out.print(_intermediateResultFile.getAbsolutePath());
        //searcher.setIntermediateFile(_intermediateResultFile.getAbsolutePath());

        Func1<Network, Double> scorer = getScoreFunction(gtsForNetworkInference, species2alleles, gtCorrespondence, singleAlleleSpecies);
        searcher.search(speciesNetwork, scorer, 1, _numRuns, _maxExaminations, _maxFailure, _optimizeBL, resultList); // search starts here



        // Move the output here. Assume that there are enough rounds so that every network
        // is populated with a real one until the maxReticulations (real reticulation node number + 2)

        // Find numPopulated, which is normally maxReticulations. However, in rare situations,
        // it can be as less than maxReticulations.



        int numPopulated = -1;
        for (int j=0; j<=9; j++) {
            if (_optimalScoreswithReticulations[j] == Double.NEGATIVE_INFINITY) {
                numPopulated = j-1;
                break;
            }
        }
        if(_printDetails) {
            System.out.println("numPopulated = " + numPopulated);
            for (int j = 0; j <= numPopulated; j++) {
                System.out.println("_optimalNetworkswithReticulations[" + j + "] =" + _optimalNetworkswithReticulations[j].toString());
                System.out.println("_optimalScoreswithReticulations[" + j + "] =" + _optimalScoreswithReticulations[j]);
            }
        }
        int correctNumReticulations = doCrossValidation(originalGTs, species2alleles, gtsForNetworkInference, gtCorrespondence, maxReticulations, singleAlleleSpecies);

        resultList.clear();
        resultList.add(new Tuple<Network, Double>(_optimalNetworkswithReticulations[correctNumReticulations], _optimalScoreswithReticulations[correctNumReticulations]));
    }

    private int doCrossValidation(List<List<MutableTuple<Tree,Double>>> originalGTs, Map<String,List<String>> species2alleles, List<Tree> distinctBinaryTrees, List<Tuple<MutableTuple<Tree, Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList, int maxReticulations, Set<String> singleAlleleSpecies){
               /* Cross validation
            1. The external while loop is for the number of reticulations.  Undetermined, go up until stoppage criteria are
               met, or exceeding the maximum reticulations limit.

            2. The internal for loop is for cross validation. K rounds for K-fold cross validation using Luay's new
                algorithm.

            3. myContainer_t and myContainer_v hold the tsublist (training sublist) and vsublist (validation sublist)
               MSE results for every sublist and every reticulation node. myContainer2_t and myContainer2_v hold
               the total training MSE and total validation MSE results for every reticulation node.

            4: The process of the internal loop.
            a. Given certain k, retrieve the model, train it to optimal branch lengths and probabilities with the
               tsublist. Return ML value, probList
               (for each distinct gene tree in the all gene tree list),and the optimized network.
            b. Use probList to compute training MSE and save it into myContainer_t.
            c. Use probList to compute validation MSE and save it into myContainer_r.

            5. If the weight of a tree is zero in a sublist, do not pass tuple information
                in nbTreeAndCountAndBinaryIDListForCV to nbTreeAndCountAndBinaryIDList of that sublist.

            6. The external process outside the double loop.
               MSE_vsublist = \sum{vsublist error squared}/numDistinctTrees in vsublist.
               MSE_tsublist = \sum{tsublist error squared}/numDistinctTrees in tsublist;
         */

        List<Tree> gts = new ArrayList<Tree>();
        for(List<MutableTuple<Tree,Double>> list: originalGTs){
            for(MutableTuple<Tree,Double> tuple: list) {
                Tree tr = tuple.Item1;
                ((STINode<Double>) tr.getRoot()).setData(tuple.Item2);
                gts.add(tr);
            }
        }

        List<Tuple3<Tree, List<Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDListForCV =
                new ArrayList<Tuple3<Tree, List<Double>, Set<Integer>>>();  // New data structure
        summarizeGeneTreestoKfolds(gts, distinctBinaryTrees, nbTreeAndCountAndBinaryIDListForCV);

        int ret = 0; // initial and current value of reticulation nodes
        int correctK = -1; // the correct number of reticulation nodes

        List<Double> myContainer_t = new ArrayList<Double>();  // container for tsublists (1 to K) MSE for k=i
        List<Double> myContainer_v = new ArrayList<Double>();  // container for vsublists (1 to K) MSE for k=i
        List<Double> myContainer2_t = new ArrayList<Double>(); // container for total MSE for all reticulations
        List<Double> myContainer2_v = new ArrayList<Double>(); // container for total MSE for all reticulations

        while (ret <= maxReticulations){   // real 2, consider 0 through 4 to pinpoint 3 if necessary
            if(_printDetails){
                System.out.println(); // restart from a line
            }

            // Loop for K tsublists to train the network models
            for (int s=_numFolds+1; s<= _numFolds*2; s++) {
                if(_printDetails){
                    System.out.println("Current reticulation = " + ret + "  s = "+s);
                }

                // Retrieve the network model from _optimalNetworkwithReticulations with the right reticulation
                Network network = _optimalNetworkswithReticulations[ret];

                // Transfer info from nbTreeAndCountAndBinaryIDListForCV to a new
                // nbTreeAndCountAndBinaryIDList used for this tsublist
                nbTreeAndCountAndBinaryIDList.clear(); // reuse the old one but with empty array
                for (Tuple3<Tree,List<Double>,Set<Integer>> tuple3: nbTreeAndCountAndBinaryIDListForCV) {
                    if (tuple3.Item2.get(s) > 0) { // weight > 0 --- zero weight terms not included
                        Tuple<MutableTuple<Tree,Double>, Set<Integer>> tuple =
                                new Tuple(new MutableTuple(tuple3.Item1, tuple3.Item2.get(s)),tuple3.Item3);
                        nbTreeAndCountAndBinaryIDList.add(tuple);
                    }
                }

                // Use the correct part of the gts for the model.
                int numVsublistTrees = gts.size()/_numFolds;
                List<Tree> gts3 = new ArrayList<Tree>();
                int s1 = s-_numFolds; // the section of gts not used --- _numFolds*(s1-1) to _numFolds*s1 (exclusive)
                if (s1 == 1) {
                    gts3.addAll(gts.subList(numVsublistTrees,gts.size()));
                }
                else {
                    gts3.addAll(gts.subList(0,numVsublistTrees*(s1-1)));
                    gts3.addAll(gts.subList(numVsublistTrees*s1,gts.size()));
                }

                // Given the network, use Brent to optimize its lengths and probabilities.
                // It modifies the species network to an optimal one during the Brent process.
                double score = computeLikelihood(network, species2alleles, distinctBinaryTrees, nbTreeAndCountAndBinaryIDList, singleAlleleSpecies);

                // Compute theoretical probabilities and real frequencies with trained network and validation data set.

                // Get theoretical probability/frequency of each distinct gene tree from the tsublist trained network
                GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
                double[] probList = new double[distinctBinaryTrees.size()];
                // The computed probability list is for the distinct binary trees in the whole list.
                gtp.calculateGTDistribution(network, distinctBinaryTrees, species2alleles, probList);

                // Compute real frequency for each distinct binary tree (distinctTrees) in tsublist and vsublist respectively.
                //
                // The total weight of a non-binary tree, its binary refinements,
                // and the theoretical probability of each distinct binary tree are known.
                // Pick the binary trees with highest theoretical probability and put
                // the averaged weight on them.

                // 1. Work on nbTreeAndCountAndBinaryIDList (tsublist).
                List<Double> realWeightsDistinctTrees_tsublist = new ArrayList<Double>();
                for (int i=0; i<=distinctBinaryTrees.size()-1;i++)
                    realWeightsDistinctTrees_tsublist.add(0.0); // initialize weights to 0.0

                for (Tuple<MutableTuple<Tree,Double>,Set<Integer>> triple: nbTreeAndCountAndBinaryIDList) {
                    // For each non-binary tree, find the set of binary gene trees
                    // with the highest probability values.

                    // Get maxProb value
                    double maxProb = 0.0;
                    for(int id: triple.Item2) // id is the index in discreteTrees and in probList
                        maxProb = Math.max(maxProb, probList[id]);  // get maxProb value

                    // Get the number of corresponding binary trees in distinctTrees that have
                    // maxProb and their id's in discreteTrees.
                    int maxProbCount = 0;
                    List<Integer> maxProbBinaryIDList = new ArrayList<Integer>();
                    for (int id: triple.Item2) {
                        if (Math.abs(maxProb - probList[id])/maxProb<1e-10) {
                            // Essentially maxProb = probList.get(id)
                            maxProbCount++;            // number of trees have maxProb
                            maxProbBinaryIDList.add(id);   // the tree ids that have maxProb
                        }
                    }

                    // Compute the average weight.
                    double weight = triple.Item1.Item2;  // get the non-binary tree weight
                    double averageWeight = weight / (double) maxProbCount;  // average weight to the GTs.

                    // Add weight to realWeightsDistinctTrees_tsublist
                    for (int id: maxProbBinaryIDList) {
                        double origWeight = realWeightsDistinctTrees_tsublist.get(id);
                        realWeightsDistinctTrees_tsublist.set(id,origWeight+averageWeight);
                    }
                }  // for each non-binary tree

                // Convert real weights on distinct trees to their frequencies in tsublist
                // Find the sum of realWeightsDistinctTrees_tsublist and the number of nonzero terms -- count
                double sum = 0;
                int count = 0;   // how many distinct trees have representations in tsublist.
                for (int i=0; i<=realWeightsDistinctTrees_tsublist.size()-1;i++) {
                    if (realWeightsDistinctTrees_tsublist.get(i)>0.0) {  // exclude distinct trees with zero weights
                        count++;  // nonzero terms
                        sum += realWeightsDistinctTrees_tsublist.get(i);
                    }
                }

                // Generate the realFreqsDistinctTrees_tsublist
                List<Double> realFreqsDistinctTrees_tsublist = new ArrayList<Double>();
                for (int i=0; i<=realWeightsDistinctTrees_tsublist.size()-1;i++)
                    realFreqsDistinctTrees_tsublist.add(realWeightsDistinctTrees_tsublist.get(i)/sum);

                // Compute the MSE for the tsublist and store MSE into myContainer_t
                double SE = 0;
                double oneError;
                for (int i=0; i<=realFreqsDistinctTrees_tsublist.size()-1;i++) {
                    if (realFreqsDistinctTrees_tsublist.get(i)>0.0) { // consider only types existent in tsublist
                        oneError = realFreqsDistinctTrees_tsublist.get(i)-probList[i];
                        SE += oneError*oneError;
                    }
                }
                double MSE = SE/count;
                myContainer_t.add(MSE);


                // 2. Work on nbTreeAndCountAndBinaryIDList_vsublist.
                // Transfer info from nbTreeAndCountAndBinaryIDListForCV
                List<Tuple3<Tree,Double, Set<Integer>>> nbTreeAndCountAndBinaryIDList_vsubList = new ArrayList<Tuple3<Tree,Double, Set<Integer>>>();

                for (Tuple3<Tree,List<Double>,Set<Integer>> tuple: nbTreeAndCountAndBinaryIDListForCV) {
                    if (tuple.Item2.get(s-_numFolds) > 0) { // weight > 0, only those present in vsublist
                        Tuple3<Tree,Double,Set<Integer>> triple =
                                new Tuple3(tuple.Item1,
                                        tuple.Item2.get(s-_numFolds),
                                        tuple.Item3);
                        nbTreeAndCountAndBinaryIDList_vsubList.add(triple);
                    }
                }

                List<Double> realWeightsDistinctTrees_vsublist = new ArrayList<Double>();
                for (int i = 0; i <= distinctBinaryTrees.size()-1;i++)
                    realWeightsDistinctTrees_vsublist.add(0.0); // initialize weights to 0.0


                for (Tuple3<Tree,Double,Set<Integer>> triple: nbTreeAndCountAndBinaryIDList_vsubList) {
                    // For each non-binary tree, find the set of binary gene trees
                    // with the highest probability values.

                    // Get maxProb value
                    double maxProb = 0.0;
                    for (int id: triple.Item3){ // id is the index in discreteTrees and in probList
                        maxProb = Math.max(maxProb, probList[id]); // get maxProb value
                    }

                    // Get the number of corresponding binary trees in distinctTrees that have
                    // maxProb and their id's in discreteTrees.
                    int maxProbCount = 0;
                    List<Integer> maxProbBinaryIDList = new ArrayList<Integer>();
                    for (int id: triple.Item3) {
                        if (Math.abs(maxProb - probList[id])/maxProb<1e-10) {
                            // Essentially maxProb = probList.get(id)
                            maxProbCount++;            // number of trees have maxProb
                            maxProbBinaryIDList.add(id);   // the tree ids that have maxProb
                        }
                    }

                    // Compute the average weight.
                    double weight = triple.Item2;  // get the non-binary tree weight
                    double averageWeight = weight / (double) maxProbCount;  // average weight to the GTs.

                    // Add weight to realWeightsDistrinctTrees_vsublist
                    for (int id: maxProbBinaryIDList) {
                        double origWeight = realWeightsDistinctTrees_vsublist.get(id);
                        realWeightsDistinctTrees_vsublist.set(id,origWeight+averageWeight);
                    }
                }

                // Convert real weights on distinct trees to their frequencies in vsublist
                // Find the sum of realWeightsDistinctTrees_vsublist
                sum = 0;
                count = 0;   // how many distinct trees have representations in vsublist.
                for (int i=0; i<=realWeightsDistinctTrees_vsublist.size()-1;i++) {
                    if (realWeightsDistinctTrees_vsublist.get(i)>0.0) {
                        count++;
                        sum += realWeightsDistinctTrees_vsublist.get(i);
                    }
                }

                List<Double> realFreqsDistinctTrees_vsublist = new ArrayList<Double>();
                for (int i=0; i<=realWeightsDistinctTrees_vsublist.size()-1;i++)
                    realFreqsDistinctTrees_vsublist.add(realWeightsDistinctTrees_vsublist.get(i)/sum);

                // Compute the MSE for the vsublist and store MSE into myContainer_v
                SE = 0;
                for (int i=0; i<=realFreqsDistinctTrees_vsublist.size()-1;i++)
                    if (realFreqsDistinctTrees_vsublist.get(i)>0.0) {
                        oneError = realFreqsDistinctTrees_vsublist.get(i)-probList[i];
                        SE += oneError*oneError;
                    }
                MSE = SE/count;
                myContainer_v.add(MSE);

                // Note: Since we are only using a sub-list, it may not have
                // all the distinct gene trees as in the distinctTrees. Therefore,
                // the real frequency should be with respect to these trees that
                // appeared in the sublist.
            }  // for, K-fold rotation

            // Compute the total MSE corresponding to the current k value
            double SE_t = 0;
            double SE_v = 0;
            for (double entry: myContainer_t)
                SE_t += entry;
            myContainer_t.clear();

            for (double entry: myContainer_v)
                SE_v += entry;
            myContainer_v.clear();

            myContainer2_t.add(SE_t/_numFolds);  // store MSE_t
            myContainer2_v.add(SE_v/_numFolds);  // store MSE_v

            if(_printDetails) {
                for (double entry : myContainer2_t) {
                    System.out.println("In myContainer2_t:" + entry);
                }
                for (double entry : myContainer2_v) {
                    System.out.println("In myContainer2_v:" + entry);
                }
            }

            // Judgement criteria: going up.
            if (myContainer2_t.size() >=2) {

                // If MSE_t not always going down, something wrong. Stop.
                for (int i=0; i<=myContainer2_t.size()-2;i++) {
                    double MSE1_t = myContainer2_t.get(i);
                    double MSE2_t = myContainer2_t.get(i+1);
                    if (MSE2_t >= MSE1_t) {
                        System.out.println("Warning: MSE increases in the training set.");
                    }
                }

                // Found up point.
                for (int i=0; i<=myContainer2_v.size()-2; i++) {
                    double MSE1_v = myContainer2_v.get(i);
                    double MSE2_v = myContainer2_v.get(i+1);
                    if (MSE2_v > MSE1_v) {
                        System.out.println("Number of reticulation nodes = " + i);
                        correctK = i;
                        return correctK;
                    }
                }
            }
            ret++;
        } // while

        // Found slow declining point for the whole curve
        for (int i=0; i<=myContainer2_v.size()-3; i++) {
            double MSE1_v = myContainer2_v.get(i);
            double MSE2_v = myContainer2_v.get(i+1);
            double MSE3_v = myContainer2_v.get(i+2);
            if (MSE1_v > MSE2_v && MSE2_v > MSE3_v &&
                    (MSE1_v-MSE2_v) > (MSE2_v-MSE3_v) &&
                    (MSE2_v-MSE3_v)/MSE3_v< 3e-2) { // turning point found by sticky area criteria
                if(_printDetails){
                    System.out.println("Turning point found within 3% MSE after larger drop in MSE. The right reticulation number is " + (i+1));
                }
                correctK = i+1;
                return correctK;

            }
        }

        // Now I have nothing captured.
        if(_printDetails){
            System.out.println("Tested 0, 1, ..., and maxReticulations. No turning point found.");
        }

        return -1;

    }


    private void summarizeGeneTreestoKfolds(List<Tree> originalGTs,
                                            List<Tree> distinctGTs,
                                            List<Tuple3<Tree, List<Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDListForCV){

        // Input:
        // 1. originalGTs: input gene trees
        // 2. distinctGTs: distinct binary gene trees in the originalGTs set
        // 3. nbTreeAndCountAndBinaryIDListForCV
        // Output:
        //  1. originalGTs: same
        //  2. distinctGTs: established.
        //  3. nbTreeAndCountAndBinaryIDListFor CV: established.
        // Note that I changed the definition of nbTreeAndCountAndBinaryIDListForCV.

        // Make sure the total number of gene trees numGT mod K = 0. Otherwise throw an exception.
        int numOriginalGTs = originalGTs.size();
        try {
            if (numOriginalGTs % _numFolds != 0)
                throw new RuntimeException("The number of gene tree is not an integer multiple of K folds");
        }
        catch (RuntimeException e) {
            System.out.println("In order for cross validation to work, you need to provide the number of " +
                    "gene trees in an integer multiple of the number of folds, whose default is 10.");
            System.exit(-1);
        }

        // Start the 1st major loop
        // Check types of gene trees and find the number of each type in originalGTs
        for(Tree tr: originalGTs){

            // The structure of nbTreeAndCountAndBinaryIDList. A triple = (nbTree, CountList, BinaryIDList)
            //
            // nbTree is a gene tree, it could be non-binary.  In the original program,
            //
            // Count = the number of effective occurrences (all weights added together) of the nbTree in
            // originalGTs.
            //
            // Expand Count to CountList in nbTreeAndCountAndBinaryIDListForCV into an array of 2K+1 integers.
            //
            // The 1st double = effective occurrence count in the originalOriginalGTs.
            // The 2nd double = effective occurrence count in the 1st vsubList
            // ...
            // The (K+1)th double = effective occurrence count in the (K+1)th vsublist.
            // The (K+1+1)th double = effective occurrence count in the 1st tsubList
            // The (K+1+2)th double = effective occurrence count in the 2nd tsubList
            // ...
            // The (2K+1)th double = effective occurrence count in the Kth tsubList
            //
            // Note:
            // 1. In this loop, OriginalGTs contains all trees.
            // 2. Since the ith double + the (i+K)th double = the 1st double in CountList,
            //    use the 1st double - the ith double to get the (i+K)th double
            //
            // BinaryIDList = a list of the binaryIDs -- the nbTree's binary resolution trees' indices in
            // distinctTrees.
            //
            // For each (non-binary) gene tree, we have its effective occurrences in the whole list,
            // in each sub-list, and the indices of its refined gene trees in distinctGTs.

            Double weight = ((STINode<Double>)tr.getRoot()).getData();
            if (weight == null) {
                weight=1.0;
            }
            int index = 0;     // the index of nbTreeAndCountAndBinaryIDListForCV
            Tuple3<Tree, List<Double>, Set<Integer>> newTuple = null;  // initial tuple

            for(Tuple3<Tree, List<Double>, Set<Integer>> triple: nbTreeAndCountAndBinaryIDListForCV){
                // Loop for all tuples in nbTreeAndCountAndBinaryIDListForCV up to now.

                if(Trees.haveSameRootedTopology(tr, triple.Item1)){
                    // The same topology tree, tr matches a triple's tree in nbTreeAndCountAndBinaryIDListForCV.
                    //
                    // Update the triple entry and give it to newTuple, and jump out of the for loop
                    //
                    triple.Item2.set(0, triple.Item2.get(0)+weight);  // update the total occurrence stored in Item2
                    newTuple = new Tuple3<Tree, List<Double>, Set<Integer>>(tr, triple.Item2, triple.Item3);
                    // Create a newTuple = an updated version of the one saved in
                    // nbTreeAndCountAndBinaryIDListForCV

                    break;  // when breaking out, the index points to this item
                }
                index++;  // if no match, increment index to be the next item.

            } // for loop for all tuples in nbTreeAndCountAndBinaryIDListForCV up to now.

            if(newTuple!=null){ // tr Matched triple and jumped out case,
                // put updated newTuple into nbTreeAndCountAndBinaryIDListForCV at index position
                // that is, to update nbTreeAndCountAndBinaryIDListForCV

                nbTreeAndCountAndBinaryIDListForCV.set(index, newTuple);
            }
            else{  // Not match case, tr did not match a triple's Item 1 in nbTreeAndCountAndBinaryIDListForCV.
                // Create a new binaryID for each of tr's binary resolutions and
                // put them into binaryIDs.

                Set<Integer> binaryIDs = new HashSet<Integer>();
                for(Tree btr: Trees.getAllBinaryResolution(tr)){
                    // loop for all binary resolution of tr, tr could be non-binary
                    // btr is one of its binary resolution trees.
                    // distinctGTs consists of all binary trees.

                    index = 0;  // this index represents the index in distinctGTs (technique)
                    boolean exist = false;
                    for(Tree exTr: distinctGTs){   // loop for all elements in distinctGTs
                        if(Trees.haveSameRootedTopology(btr, exTr)){
                            // a tr's binary resolution matches a distinctGT
                            binaryIDs.add(index); // add this distinctGT's index to binaryID
                            exist = true;   // This tree is contained in distinctGTs already
                            break;  // jump out of the loop for all elements in distinctGTs
                        }
                        index++; // no match, increment index to next one in distinctGTs
                        // if the loop ends normally without break, index = 1 larger
                        // than the index of the last element in distinctGTs

                    }  // loop for all elements in distinctGTs

                    if(!exist){ // no match is found, this binary tree is not in distinctGTs
                        distinctGTs.add(btr);  // add this binary tree into distinctGTs
                        binaryIDs.add(index);
                        // add the new element's
                        // index in the new distinctGT to binaryIDs
                    }
                }   // loop for all binary resolution of tr

                // Create a new Tuple to be put into nbTreeAndCountAndBinaryIDListForCV
                // Create the Countlist
                ArrayList<Double> temp = new ArrayList<Double>(_numFolds*2+1);
                temp.add(0,weight);  // the first element of the CountList
                for (int i=1; i<=_numFolds*2;i++)  // the other elements of the CountList
                    temp.add(0.0);

                newTuple = new Tuple3<Tree, List<Double>, Set<Integer>>(tr, temp, binaryIDs);

                nbTreeAndCountAndBinaryIDListForCV.add(newTuple);
                // This new tr is a first encounter. Together with its binaryID's in
                // distinctGTs, a newTuple is created and is put into nbTreeAndCountAndBinaryIDListForCV
            }  // else for no match case
        } // loop for originalGTs


        // Start the 2nd major loop
        //
        // Divide all gene trees into K sub-lists.
        // Find the number of each type of gene trees in each sub-list of originalGTs
        //
        // For each sub-list, compute gene tree numbers of each type.
        // For each complementary sub-list, compute gene tree numbers of each type.
        //
        // Loop 1: from 1 to K and create a sub-list of originalGT each time.
        // Loop 2: for that sub-list (and the corresponding complementary sub-list) created in loop 1.

        int numASubListOriginalGTs = numOriginalGTs/_numFolds;

        // Loop for all vsublists, compute the correct frequencies of each tree
        for (int i=0; i<=_numFolds-1;i++) {
            List<Tree> aSubListOriginalGTs = originalGTs.subList(numASubListOriginalGTs*i,
                    numASubListOriginalGTs*(i+1));  // create vsublist

            for(Tree tr: aSubListOriginalGTs){ // loop for the sub-list of original gene trees.
                Double weight = ((STINode<Double>)tr.getRoot()).getData();
                if (weight == null) {
                    weight = 1.0;
                }
                int index = 0;     // the index of nbTreeAndCountAndBinaryIDListForCV
                Tuple3<Tree, List<Double>, Set<Integer>> newTuple = null;  // initial tuple

                for(Tuple3<Tree, List<Double>, Set<Integer>> triple: nbTreeAndCountAndBinaryIDListForCV){
                    // Loop for all tuples in nbTreeAndCountAndBinaryIDListForCV.

                    if(Trees.haveSameRootedTopology(tr, triple.Item1)){
                        // The same topology tree, tr matches a triple's tree in nbTreeAndCountAndBinaryIDList.
                        // Update CountList and give it to newTuple, and jump out of this loop

                        double temp1 = triple.Item2.get(i+1);  // Get the original value of the i+1 th term
                        // in CountList
                        temp1 +=weight;  // increment temp1
                        triple.Item2.set(i+1,temp1); // increment the (i+1)th element in CountList

                        newTuple = new Tuple3<Tree, List<Double>, Set<Integer>>(tr, triple.Item2, triple.Item3);
                        // create a newTuple = an updated version of the one saved in nbTreeAndCountAndBinaryIDListForCV

                        break;  // when breaking out, the index points to this item
                    }
                    index++;  // if no match, increment index to be the next item.

                } // for loop for all tuples in nbTreeAndCountAndBinaryIDListForCV.
                // Note: In this round every tree must find its counterpart in
                // nbTreeAndCountAndBinaryIDListForCV. Otherwise throw an exception.

                if(newTuple!=null){ // tr Matched triple and jumped out case,
                    // put updated newTuple into nbTreeAndCountAndBinaryIDList at index position
                    // that is, to update nbTreeAndCountAndBinaryIDList

                    nbTreeAndCountAndBinaryIDListForCV.set(index, newTuple);
                }
                else  // Not match case, throw error
                    throw new RuntimeException("New gene tree topology is found. Error!");
            } // loop for the vsublist of original gene trees.

        }  // Loop for all vsublists, compute the correct frequencies

        // Update the tsublist related frequences
        for(Tuple3<Tree, List<Double>, Set<Integer>> triple: nbTreeAndCountAndBinaryIDListForCV){
            double totalWeight = triple.Item2.get(0);  // Get the original value of the 1st term in CountList
            for (int i=1; i<=_numFolds;i++) {
                double vsublistWeight = triple.Item2.get(i);
                triple.Item2.set(i+_numFolds,totalWeight-vsublistWeight); // Update the (i+1+K)th element in CountList
            }
        }
    }


}
