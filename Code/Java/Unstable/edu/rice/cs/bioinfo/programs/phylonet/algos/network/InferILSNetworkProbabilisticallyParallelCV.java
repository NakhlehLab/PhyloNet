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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.RichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.phylogenetics.FindRoot;
import edu.rice.cs.bioinfo.library.phylogenetics.GetDirectSuccessors;
import edu.rice.cs.bioinfo.library.phylogenetics.GetInDegree;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkNeighbourhoodRandomWalkGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours.AllNeighboursHillClimberFirstBetter;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.io.ByteArrayInputStream;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: jianrongdong
 * Date: 7/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferILSNetworkProbabilisticallyParallelCV extends MDCOnNetworkYFFromRichNewickJung {
    protected Network[] _optimalNetworkswithReticulations;
    protected double[] _optimalScoreswithReticulations;
    protected int _numMultipleRuns; // the number of multiple runs, either with same starting point or not
    protected int _maxRounds;
    protected int _maxTryPerBranch;
    protected double _improvementThreshold;
    protected double _maxBranchLength;
    protected double _Brent1;
    protected double _Brent2;
    protected Long _maxFailure;
    protected Long _maxExaminations;
    protected int _diameterLimit;
    protected Network<Object> _startNetwork;
    protected Set<String> _fixedHybrid;
    protected int _numThread;
    protected double[] _operationWeight;
    protected int _numFolds;
    protected Long _seed;
    protected boolean _printDetails = false;



    public InferILSNetworkProbabilisticallyParallelCV(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }

    public void setSearchParameter(int maxRounds,
                                   int maxTryPerBranch,
                                   double improvementThreshold,
                                   double maxBranchLength,
                                   double Brent1,
                                   double Brent2,
                                   Long maxExaminations,
                                   Long maxFailure,
                                   int diameterLimit,
                                   int parallel,
                                   Network startNetwork,
                                   Set<String> fixedHybrid,
                                   int numMultipleRuns,
                                   int numFolds,
                                   double[] operationWeight,
                                   Long seed) {
        _maxRounds = maxRounds;
        _maxTryPerBranch = maxTryPerBranch;
        _improvementThreshold = improvementThreshold;
        _maxBranchLength = maxBranchLength;
        //_maxBranchLength = 12;
        _Brent1 = Brent1;
        _Brent2 = Brent2;
        _maxExaminations = maxExaminations;
        _diameterLimit = diameterLimit;
        _startNetwork = startNetwork;
        _maxFailure = maxFailure;
        _numThread = parallel;
        _fixedHybrid = fixedHybrid;
        _numMultipleRuns = numMultipleRuns;
        _numFolds = numFolds;
        _operationWeight = operationWeight;
        _seed = seed;


    }

    private void checkNetworkWithHybrids(Network<Object> startNetwork){
        if(startNetwork == null || _fixedHybrid.size() == 0){
            return;
        }

        if(startNetwork != null){
            for(NetNode<Object> node: startNetwork.getNetworkNodes()){
                NetNode hybridSpecies = node.getChildren().iterator().next();
                if(!(hybridSpecies.isLeaf() && _fixedHybrid.contains(hybridSpecies.getName()))){
                    throw new IllegalArgumentException("The starting network contains hybrid that is not in the specified hybrid set.");
                }
            }
        }

    }

    public Tuple3<Network,Double,Integer> CV(List<MutableTuple<Tree,Double>> originalGTs, Map<String,List<String>> species2alleles,
                  int maxReticulations, int hasTried){
        _optimalNetworkswithReticulations = new Network[10];   // hold the best network so far with k reticulation nodes
        _optimalScoreswithReticulations = new double[10];     // hold the ML score with the corresponding k reticulation nodes
        int numPopulated = -1;

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

        List<MutableTuple<Tree,Double>> distinctNonBinaryTrees = new ArrayList<MutableTuple<Tree, Double>>();
        summarizeGeneTrees1(originalGTs, allele2species, distinctNonBinaryTrees);

        String startingNetwork = getStartNetwork(distinctNonBinaryTrees, species2alleles, _fixedHybrid, _startNetwork);
        //System.out.println(startingNetwork);
        List<Tree> distinctBinaryTrees = new ArrayList<Tree>();
        List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList = new ArrayList<Tuple3<Tree, Double, List<Integer>>>();
        summarizeGeneTrees2(distinctNonBinaryTrees, allele2species, distinctBinaryTrees, nbTreeAndCountAndBinaryIDList);

        List<Tree> gts = new ArrayList<Tree>();
        for(MutableTuple<Tree,Double> tuple: originalGTs){
            Tree tr = tuple.Item1;
            ((STINode<Double>)tr.getRoot()).setData(tuple.Item2);
            gts.add(tr);
        }

        // Multiple runs. During each run, the best network found with each reticulation node is compared
        // with that network inside _optimalNetworkwithReticulations. Hence we save the best networks from multiple runs.
        if(_printDetails){
            System.out.println();
        }
        for (int i=0; i< _numMultipleRuns; i++) {
            if(_printDetails){
                System.out.println("i = " + i);
            }

           // _startNetwork = null;
            DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = makeNetwork(startingNetwork);

            NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(_operationWeight, makeNode, makeEdge, _seed);
            AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);

            Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> scorer = getScoreFunction(distinctBinaryTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
            Comparator<Double> comparator = getDoubleScoreComparator();
            HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations, _maxFailure, _diameterLimit, hasTried); // search starts here
            System.gc();

            // Move the output here. Assume that there are enough rounds so that every network
            // is populated with a real one until the maxReticulations (real reticulation node number + 2)

            // Find numPopulated, which is normally maxReticulations. However, in rare situations,
            // it can be as less than maxReticulations.

            for (int j=0; j<=9; j++) {
                if (_optimalScoreswithReticulations[j] == Double.NEGATIVE_INFINITY) {
                    numPopulated = j-1;
                    break;
                }
            }
            if(_printDetails) {
                System.out.println("numPopulated = " + numPopulated);
                for (int j = 0; j <= numPopulated; j++) {
                    System.out.println("_optimalNetworkswithReticulations[" + j + "] =" + network2String(_optimalNetworkswithReticulations[j]));
                    System.out.println("_optimalScoreswithReticulations[" + j + "] =" + _optimalScoreswithReticulations[j]);
                }
            }
        }

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

        List<Tree> gts2 = new ArrayList<Tree>(gts);   // get a deep copy of gts
        List<Tuple3<Tree, List<Double>, List<Integer>>> nbTreeAndCountAndBinaryIDListForCV =
                new ArrayList<Tuple3<Tree, List<Double>, List<Integer>>>();  // New data structure
        summarizeGeneTreestoKfolds(gts2, distinctBinaryTrees, nbTreeAndCountAndBinaryIDListForCV);

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
                for (Tuple3<Tree,List<Double>,List<Integer>> tuple: nbTreeAndCountAndBinaryIDListForCV) {
                    if (tuple.Item2.get(s) > 0) { // weight > 0 --- zero weight terms not included
                        Tuple3<Tree,Double,List<Integer>> triple =
                                new Tuple3<Tree,Double,List<Integer>>(tuple.Item1, tuple.Item2.get(s),tuple.Item3);
                        nbTreeAndCountAndBinaryIDList.add(triple);
                    }
                }

                // Use the correct part of the gts for the model.
                int numVsublistTrees = gts2.size()/_numFolds;
                List<Tree> gts3 = new ArrayList<Tree>();
                int s1 = s-_numFolds; // the section of gts not used --- _numFolds*(s1-1) to _numFolds*s1 (exclusive)
                if (s1 == 1) {
                    gts3.addAll(gts2.subList(numVsublistTrees,gts2.size()));
                }
                else {
                    gts3.addAll(gts2.subList(0,numVsublistTrees*(s1-1)));
                    gts3.addAll(gts2.subList(numVsublistTrees*s1,gts2.size()));
                }

                // Given the network, use Brent to optimize its lengths and probabilities.
                // It modifies the species network to an optimal one during the Brent process.
                double score = findNonUltrametricOptimalBranchLength(network,
                                                                     distinctBinaryTrees, species2alleles,
                                                                     nbTreeAndCountAndBinaryIDList);

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

                for (Tuple3<Tree,Double,List<Integer>> triple: nbTreeAndCountAndBinaryIDList) {
                    // For each non-binary tree, find the set of binary gene trees
                    // with the highest probability values.

                    // Get maxProb value
                    double maxProb = 0.0;
                    for(int id: triple.Item3) // id is the index in discreteTrees and in probList
                        maxProb = Math.max(maxProb, probList[id]);  // get maxProb value

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
                List<Tuple3<Tree,Double, List<Integer>>> nbTreeAndCountAndBinaryIDList_vsubList = new ArrayList<Tuple3<Tree,Double, List<Integer>>>();

                for (Tuple3<Tree,List<Double>,List<Integer>> tuple: nbTreeAndCountAndBinaryIDListForCV) {
                    if (tuple.Item2.get(s-_numFolds) > 0) { // weight > 0, only those present in vsublist
                        Tuple3<Tree,Double,List<Integer>> triple =
                                new Tuple3<Tree,Double,List<Integer>>(tuple.Item1,
                                        tuple.Item2.get(s-_numFolds),
                                        tuple.Item3);
                        nbTreeAndCountAndBinaryIDList_vsubList.add(triple);
                    }
                }

                List<Double> realWeightsDistinctTrees_vsublist = new ArrayList<Double>();
                for (int i = 0; i <= distinctBinaryTrees.size()-1;i++)
                    realWeightsDistinctTrees_vsublist.add(0.0); // initialize weights to 0.0


                for (Tuple3<Tree,Double,List<Integer>> triple: nbTreeAndCountAndBinaryIDList_vsubList) {
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
                        return new Tuple3<Network, Double, Integer>(_optimalNetworkswithReticulations[correctK], _optimalScoreswithReticulations[correctK], correctK);
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
                return new Tuple3<Network, Double, Integer>(_optimalNetworkswithReticulations[correctK], _optimalScoreswithReticulations[correctK], correctK);

            }
        }

        // Now I have nothing captured.
        if(_printDetails){
            System.out.println("Tested 0, 1, ..., and maxReticulations. No turning point found.");
        }
        return new Tuple3<Network, Double, Integer>(_optimalNetworkswithReticulations[correctK], _optimalScoreswithReticulations[correctK], correctK);
    }

    // Modified from summarizeGeneTrees
    private void summarizeGeneTreestoKfolds(List<Tree> originalGTs,
                                            List<Tree> distinctGTs,
                                            List<Tuple3<Tree, List<Double>, List<Integer>>> nbTreeAndCountAndBinaryIDListForCV){

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
            Tuple3<Tree, List<Double>, List<Integer>> newTuple = null;  // initial tuple

            for(Tuple3<Tree, List<Double>, List<Integer>> triple: nbTreeAndCountAndBinaryIDListForCV){
                // Loop for all tuples in nbTreeAndCountAndBinaryIDListForCV up to now.

                if(Trees.haveSameRootedTopology(tr, triple.Item1)){
                    // The same topology tree, tr matches a triple's tree in nbTreeAndCountAndBinaryIDListForCV.
                    //
                    // Update the triple entry and give it to newTuple, and jump out of the for loop
                    //
                    triple.Item2.set(0, triple.Item2.get(0)+weight);  // update the total occurrence stored in Item2
                    newTuple = new Tuple3<Tree, List<Double>, List<Integer>>(tr, triple.Item2, triple.Item3);
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

                List<Integer> binaryIDs = new ArrayList<Integer>();
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

                newTuple = new Tuple3<Tree, List<Double>, List<Integer>>(tr, temp, binaryIDs);

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
                Tuple3<Tree, List<Double>, List<Integer>> newTuple = null;  // initial tuple

                for(Tuple3<Tree, List<Double>, List<Integer>> triple: nbTreeAndCountAndBinaryIDListForCV){
                    // Loop for all tuples in nbTreeAndCountAndBinaryIDListForCV.

                    if(Trees.haveSameRootedTopology(tr, triple.Item1)){
                        // The same topology tree, tr matches a triple's tree in nbTreeAndCountAndBinaryIDList.
                        // Update CountList and give it to newTuple, and jump out of this loop

                        double temp1 = triple.Item2.get(i+1);  // Get the original value of the i+1 th term
                        // in CountList
                        temp1 +=weight;  // increment temp1
                        triple.Item2.set(i+1,temp1); // increment the (i+1)th element in CountList

                        newTuple = new Tuple3<Tree, List<Double>, List<Integer>>(tr, triple.Item2, triple.Item3);
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
        for(Tuple3<Tree, List<Double>, List<Integer>> triple: nbTreeAndCountAndBinaryIDListForCV){
            double totalWeight = triple.Item2.get(0);  // Get the original value of the 1st term in CountList
            for (int i=1; i<=_numFolds;i++) {
                double vsublistWeight = triple.Item2.get(i);
                triple.Item2.set(i+_numFolds,totalWeight-vsublistWeight); // Update the (i+1+K)th element in CountList
            }
        }
    }

    protected String getStartNetwork(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles, Set<String> hybridSpecies, Network<Object> startingNetwork){
        checkNetworkWithHybrids(startingNetwork);

        if(startingNetwork == null){
            Map<String,String> allele2species = null;
            if(species2alleles!=null){
                allele2species = new HashMap<String, String>();
                for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                    String species = entry.getKey();
                    for(String allele: entry.getValue()){
                        allele2species.put(allele,species);
                    }
                }
            }

            MDCInference_DP mdc = new MDCInference_DP();
            Solution sol;
            if(allele2species==null){
                sol = mdc.inferSpeciesTree(gts, false, 1, true, true, -1).get(0);
            }
            else{
                sol = mdc.inferSpeciesTree(gts, allele2species, false, 1, true, true, -1).get(0);
            }
            Tree startingTree= Trees.generateRandomBinaryResolution(sol._st);
            startingNetwork = string2Network(startingTree.toString());
        }



        for(String hybrid: hybridSpecies){
            createHybrid(startingNetwork, hybrid);
        }

        int index = 1;
        for(NetNode<Object> node: startingNetwork.dfs()){
            if(node.getName()==null || node.getName().equals("")){
                String name;
                do{
                    name = "i" + (index++);
                }while(startingNetwork.findNode(name)!=null);
                node.setName(name);
            }
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent, Double.NaN);
                node.setParentSupport(parent, Double.NaN);
                node.setParentProbability(parent, Double.NaN);
            }
        }

        String newNetwork = network2String(startingNetwork);
        return newNetwork;

    }

    protected void createHybrid(Network<Object> network, String hybrid){
        List<Tuple<NetNode,NetNode>> edgeList = new ArrayList<Tuple<NetNode,NetNode>>();
        Tuple<NetNode,NetNode> destinationEdge = null;
        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(network)){
            for(NetNode child: node.getChildren()){
                if(child.isLeaf() && child.getName().equals(hybrid)){
                    if(node.isNetworkNode()){
                        return;
                    }
                    destinationEdge = new Tuple<NetNode, NetNode>(node, child);
                }
                else{
                    edgeList.add(new Tuple<NetNode, NetNode>(node, child));
                }
            }

        }

        int numEdges = edgeList.size();
        Tuple<NetNode,NetNode> sourceEdge = edgeList.get((int)(Math.random() * numEdges));
        NetNode insertedSourceNode = new BniNetNode();
        insertedSourceNode.adoptChild(sourceEdge.Item2, NetNode.NO_DISTANCE);
        sourceEdge.Item1.removeChild(sourceEdge.Item2);
        sourceEdge.Item1.adoptChild(insertedSourceNode, NetNode.NO_DISTANCE);
        NetNode insertedDestinationNode = new BniNetNode();
        insertedDestinationNode.adoptChild(destinationEdge.Item2, NetNode.NO_DISTANCE);
        destinationEdge.Item1.removeChild(destinationEdge.Item2);
        destinationEdge.Item1.adoptChild(insertedDestinationNode, NetNode.NO_DISTANCE);
        insertedSourceNode.adoptChild(insertedDestinationNode, NetNode.NO_DISTANCE);
    }

    protected void summarizeGeneTrees1(List<MutableTuple<Tree,Double>> originalGTs, Map<String,String> allele2species, List<MutableTuple<Tree,Double>> distinctGTs){
        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        for(MutableTuple<Tree,Double> gtTuple: originalGTs){
            for(TNode node: gtTuple.Item1.getNodes()){
                node.setParentDistance(TNode.NO_DISTANCE);
            }
            String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
            MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
            if(existingTuple==null){
                existingTuple = gtTuple;
                exp2tree.put(exp, existingTuple);
            }
            else{
                existingTuple.Item2 += gtTuple.Item2;
            }
        }

        distinctGTs.addAll(exp2tree.values());
    }


    protected void summarizeGeneTrees2(List<MutableTuple<Tree,Double>> originalGTs, Map<String,String> allele2species, List<Tree> distinctGTs, List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        Map<String,Integer> exp2ID = new HashMap<String, Integer>();
        //int temp=0;
        for(MutableTuple<Tree,Double> gt: originalGTs){
            //System.out.println(temp + " :" + distinctGTs.size());
            //temp++;
            List<Integer> binaryIDs = new ArrayList<Integer>();
            for(Tree btr: Trees.getAllBinaryResolution(gt.Item1)){
                String btrExp = Trees.getLexicographicNewickString(btr, allele2species);
                Integer index = exp2ID.get(btrExp);
                if(index==null){
                    index = distinctGTs.size();
                    distinctGTs.add(btr);
                    binaryIDs.add(index);
                    exp2ID.put(btrExp, index);
                }
                else{
                    binaryIDs.add(index);
                }
            }
            nbTreeAndCountAndBinaryIDList.add(new Tuple3<Tree, Double, List<Integer>>(gt.Item1, gt.Item2, binaryIDs));

        }
    }


    protected Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }


    protected Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                //System.out.println("Start scoring ...");
                //System.out.println(network2String(network));
                //long start = System.currentTimeMillis();
                Network<Object> speciesNetwork = networkNew2Old(network);
                //double score = findUltrametricOptimalBranchLength(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
                double score = findNonUltrametricOptimalBranchLength(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);

                // A score can be lower than the scores that I currently have in the
                // optimal score list, but it is still the best score
                // maybe for a lower reticulation network. Therefore, I just
                // compare every network found with those in the _optimalNetworkwithReticulautions
                // to update it with the best network of k reticulations found so far.

                // Save the best network with k reticulations into _optimalNetworkswithReticulations

                // Find the number of network nodes in speciesNetwork.

                Iterator<NetNode<Object>> iter = speciesNetwork.getNetworkNodes().iterator();
                int indexCount = 0;
                while (iter.hasNext()) {// find the right number of reticulation nodes
                    iter.next();
                    indexCount++;
                }
                // System.out.println("indexCount =" + indexCount);

                if (score > _optimalScoreswithReticulations[indexCount]) {  // initially -inf, already assigned double value
                    _optimalNetworkswithReticulations[indexCount] = string2Network(network2String(speciesNetwork));
                    _optimalScoreswithReticulations[indexCount] = score;
                }
                // System.out.println("score = " + score);
                // for (int i=0; i<=indexCount; i++) {
                //     System.out.println("_optimalNetworkswithReticulations["+i+"] =" + network2String(_optimalNetworkswithReticulations[i]));
                //    System.out.println("_optimalScoreswithReticulations["+i+"] =" + _optimalScoreswithReticulations[i]);
                // }
                System.gc();
                return score;
            }
        };
    }


    protected double findUltrametricOptimalBranchLength(final Network<Object> speciesNetwork, final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        boolean continueRounds = true; // keep trying to improve network

        Map<NetNode<Object>, Double> node2height = new Hashtable<NetNode<Object>, Double>();
        for(NetNode<Object> parent: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            if(parent.isLeaf()){
                node2height.put(parent, 0.0);
                continue;
            }
            double height = 0;
            for(NetNode<Object> child: parent.getChildren()){
                height = Math.max(height, node2height.get(child));
            }

            height = height + 1;


            node2height.put(parent, height);
            for(NetNode<Object> child: parent.getChildren()){
                if(child.isNetworkNode()){
                    child.setParentProbability(parent,0.5);
                }

                child.setParentDistance(parent, height - node2height.get(child));
            }
        }


        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList));  // records the GTProb of the network at all times
        final Container<Map<NetNode<Object>, Double>> node2heightContainer = new Container<Map<NetNode<Object>, Double>>(node2height);


        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();
            for(final NetNode<Object> node : edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork))
            {
                if(node.isLeaf()){
                    continue;
                }

                final Container<Double> minHeight = new Container<Double>(0.0);
                final Container<Double> maxHeight = new Container<Double>(Double.MAX_VALUE);


                for(NetNode<Object> child: node.getChildren()){
                    double childHeight = node2heightContainer.getContents().get(child);
                    minHeight.setContents(Math.max(minHeight.getContents(), childHeight));

                }

                if(!node.isRoot()){
                    for(NetNode<Object> parent: node.getParents()){
                        double parentHeight = node2heightContainer.getContents().get(parent);
                        maxHeight.setContents(Math.min(maxHeight.getContents(), parentHeight));
                    }
                }
                else{
                    maxHeight.setContents(minHeight.getContents() + _maxBranchLength);
                }

                UnivariateFunction functionToOptimize = new UnivariateFunction() {

                    public double value(double suggestedHeight) {  // brent suggests a new branch length
                        double incumbentHeight = node2heightContainer.getContents().get(node);


                        for(NetNode<Object> child: node.getChildren()){
                            child.setParentDistance(node, suggestedHeight - node2heightContainer.getContents().get(child));
                        }

                        if(!node.isRoot()){
                            for(NetNode<Object> parent: node.getParents()){
                                node.setParentDistance(parent, node2heightContainer.getContents().get(parent) - suggestedHeight);
                            }
                        }

                        double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, node, true);
                        if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                        {
                            lnGtProbOfSpeciesNetwork.setContents(lnProb);
                            node2heightContainer.getContents().put(node, suggestedHeight);

                        }
                        else  // didn't improve, roll back change
                        {
                            for(NetNode<Object> child: node.getChildren()){
                                child.setParentDistance(node, incumbentHeight - node2heightContainer.getContents().get(child));
                            }
                            if(!node.isRoot()){
                                for(NetNode<Object> parent: node.getParents()){
                                    node.setParentDistance(parent, node2heightContainer.getContents().get(parent) - incumbentHeight);
                                }
                            }
                        }
                        return lnProb;
                    }
                };
                BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                try
                {
                    optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, minHeight.getContents(), maxHeight.getContents());
                }
                catch(TooManyEvaluationsException e) // _maxAssigmentAttemptsPerBranchParam exceeded
                {
                }

                computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, node, true);

            }

            for(final NetNode<Object> child : speciesNetwork.getNetworkNodes()) // find every hybrid node
            {


                Iterator<NetNode<Object>> hybridParents = child.getParents().iterator();
                final NetNode hybridParent1 = hybridParents.next();
                final NetNode hybridParent2 = hybridParents.next();

                UnivariateFunction functionToOptimize = new UnivariateFunction() {
                    public double value(double suggestedProb) {
                        double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);
                        child.setParentProbability(hybridParent1, suggestedProb);
                        child.setParentProbability(hybridParent2, 1.0 - suggestedProb);

                        double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, false);

                        if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // change improved GTProb, keep it
                        {

                            lnGtProbOfSpeciesNetwork.setContents(lnProb);
                        }
                        else // change did not improve, roll back
                        {

                            child.setParentProbability(hybridParent1, incumbentHybridProbParent1);
                            child.setParentProbability(hybridParent2, 1.0 - incumbentHybridProbParent1);
                        }
                        return lnProb;
                    }
                };
                BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                try
                {
                    optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, 0, 1.0);
                }
                catch(TooManyEvaluationsException e)  // _maxAssigmentAttemptsPerBranchParam exceeded
                {
                }
                computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, false);



            }

            if( ((double)lnGtProbOfSpeciesNetwork.getContents()) == lnGtProbLastRound)  // if no improvement was made wrt to last around, stop trying to find a better assignment
            {
                continueRounds = false;
            }
            else if (lnGtProbOfSpeciesNetwork.getContents() > lnGtProbLastRound) // improvement was made, ensure it is large enough wrt to improvement threshold to continue searching
            {

                double improvementPercentage = Math.pow(Math.E, (lnGtProbOfSpeciesNetwork.getContents() - lnGtProbLastRound)) - 1.0;  // how much did we improve over last round
                //System.out.println(improvementPercentage + " vs. " + _improvementThreshold);
                if(improvementPercentage < _improvementThreshold  )  // improved, but not enough to keep searching
                {
                    continueRounds = false;
                }
            }
            else
            {
                throw new IllegalStateException("Should never have decreased prob.");
            }
        }

        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            if(node.isLeaf()){
                if(species2alleles == null || species2alleles.get(node.getName()).size()<2){
                    node.setParentDistance(node.getParents().iterator().next(), Double.NaN) ;
                }
            }
        }

        System.out.println(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList) + " vs. " + lnGtProbOfSpeciesNetwork.getContents());
        return lnGtProbOfSpeciesNetwork.getContents();
    }



    protected double findNonUltrametricOptimalBranchLength(final Network<Object> speciesNetwork, final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        boolean continueRounds = true; // keep trying to improve network


        for(NetNode<Object> node: speciesNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent,1.0);
                if(node.isNetworkNode()){
                    //node.setParentDistance(parent,0.0);
                    node.setParentProbability(parent, 0.5);
                }
            }
        }

        //long start = System.currentTimeMillis();
        double initalProb;
        if(_numThread > 1){
            initalProb = computeProbabilityParallel(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
        }
        else{
            initalProb = computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
        }
        //System.out.print("\n"+(System.currentTimeMillis()-start));
        //computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);

        //System.out.println();
        //System.out.println(network2String(speciesNetwork));
        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(initalProb);  // records the GTProb of the network at all times
        final Container<Integer> callCount = new Container<Integer>(0);


        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
            /*
            * Prepare a random ordering of network edge examinations each of which attempts to change a branch length or hybrid prob to improve the GTProb score.
            */

            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();
            List<Proc> assigmentActions = new ArrayList<Proc>(); // store adjustment commands here.  Will execute them one by one later.




            for(final NetNode<Object> parent : edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork))
            {

                for(final NetNode<Object> child : parent.getChildren())
                {
                    if(child.isLeaf()){
                        if(species2alleles==null || species2alleles.get(child.getName()).size()<2){
                            continue;
                        }
                    }

                    assigmentActions.add(new Proc()
                    {
                        public void execute()
                        {

                            UnivariateFunction functionToOptimize = new UnivariateFunction() {
                                public double value(double suggestedBranchLength) {
                                    callCount.setContents(callCount.getContents()+1);
                                    double incumbentBranchLength = child.getParentDistance(parent);

                                    child.setParentDistance(parent, suggestedBranchLength);

                                    double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, parent);
                                    //System.out.println("Changing branch ("+parent.getName()+","+child.getName()+") to " + suggestedBranchLength);
                                    //System.out.println(network2String(speciesNetwork)+": " + lnProb);

                                    if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                                    {
                                        lnGtProbOfSpeciesNetwork.setContents(lnProb);

                                    }
                                    else  // didn't improve, roll back change
                                    {
                                        child.setParentDistance(parent, incumbentBranchLength);
                                    }
                                    return lnProb;
                                }
                            };
                            BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                            try
                            {
                                optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, Double.MIN_VALUE, _maxBranchLength);
                            }
                            catch(TooManyEvaluationsException e) // _maxAssigmentAttemptsPerBranchParam exceeded
                            {
                            }

                            computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, parent);
                        }
                    });
                }
            }


            for(final NetNode<Object> child : speciesNetwork.getNetworkNodes()) // find every hybrid node
            {

                Iterator<NetNode<Object>> hybridParents = child.getParents().iterator();
                final NetNode hybridParent1 = hybridParents.next();
                final NetNode hybridParent2 = hybridParents.next();

                assigmentActions.add(new Proc()
                {
                    public void execute()
                    {
                        UnivariateFunction functionToOptimize = new UnivariateFunction() {
                            public double value(double suggestedProb) {
                                callCount.setContents(callCount.getContents()+1);
                                double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);

                                child.setParentProbability(hybridParent1, suggestedProb);
                                child.setParentProbability(hybridParent2, 1.0 - suggestedProb);

                                double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, false);
                                //System.out.println("Changing node probability to "+ suggestedProb);
                                //System.out.println(network2String(speciesNetwork)+": " + lnProb);
                                //System.out.println(Math.abs(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList) - lnProb));
                                if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // change improved GTProb, keep it
                                {

                                    lnGtProbOfSpeciesNetwork.setContents(lnProb);
                                }
                                else // change did not improve, roll back
                                {

                                    child.setParentProbability(hybridParent1, incumbentHybridProbParent1);
                                    child.setParentProbability(hybridParent2, 1.0 - incumbentHybridProbParent1);
                                }
                                return lnProb;
                            }
                        };
                        BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                        try
                        {
                            optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, 0, 1.0);
                        }
                        catch(TooManyEvaluationsException e)  // _maxAssigmentAttemptsPerBranchParam exceeded
                        {
                        }
                        computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, false);
                        //System.out.println(network2String(speciesNetwork) + " : " + lnGtProbOfSpeciesNetwork.getContents());
                    }
                });


            }


            // add hybrid probs to hybrid edges


            for(Proc assigment : assigmentActions)   // for each change attempt, perform attempt
            {
                assigment.execute();
            }


            if( ((double)lnGtProbOfSpeciesNetwork.getContents()) == lnGtProbLastRound)  // if no improvement was made wrt to last around, stop trying to find a better assignment
            {
                continueRounds = false;
            }
            else if (lnGtProbOfSpeciesNetwork.getContents() > lnGtProbLastRound) // improvement was made, ensure it is large enough wrt to improvement threshold to continue searching
            {

                double improvementPercentage = Math.pow(Math.E, (lnGtProbOfSpeciesNetwork.getContents() - lnGtProbLastRound)) - 1.0;  // how much did we improve over last round
                if(improvementPercentage < _improvementThreshold  )  // improved, but not enough to keep searching
                {
                    continueRounds = false;
                }
            }
            else
            {
                throw new IllegalStateException("Should never have decreased prob.");
            }
        }
        //System.out.println(callCount.getContents());
        //System.out.println(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList) + " vs. " + lnGtProbOfSpeciesNetwork.getContents());
        return lnGtProbOfSpeciesNetwork.getContents();
    }



    private class MyThreadFromScratch extends Thread{
        GeneTreeProbabilityYF _gtp;
        Network _speciesNetwork;
        List<Tree> _geneTrees;
        Map<String, List<String>> _species2alleles;
        double[] _probs;


        public MyThreadFromScratch(GeneTreeProbabilityYF gtp, Network speciesNetwork, List<Tree> geneTrees, Map<String, List<String>> species2alleles, double[] probs){
            _speciesNetwork = speciesNetwork;
            _geneTrees = geneTrees;
            _species2alleles = species2alleles;
            _probs = probs;
            _gtp = gtp;
        }


        public void run() {
            _gtp.calculateGTDistribution(_speciesNetwork, _geneTrees, _species2alleles, _probs);

        }
    }


    private class MyThreadFromNonScratch extends Thread{
        Network _speciesNetwork;
        List<Tree> _gts;
        double[] _probs;
        Set<NetNode> _childNodes;
        Set<NetNode> _parentNodes;
        GeneTreeProbabilityYF _gtp;


        public MyThreadFromNonScratch(GeneTreeProbabilityYF gtp, Network speciesNetwork, List<Tree> gts, Set<NetNode> childNodes, Set<NetNode> parentNodes, double[] probs){
            _speciesNetwork = speciesNetwork;
            _gts = gts;
            _probs = probs;
            _childNodes = childNodes;
            _parentNodes = parentNodes;
            _gtp = gtp;
        }


        public void run() {
            _gtp.calculateGTDistribution(_speciesNetwork, _gts, _childNodes, _parentNodes, _probs);

        }
    }


    protected double computeProbabilityParallel(Network<Object> speciesNetwork, List<Tree> distinctTrees, Map<String, List<String>> species2alleles, List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList) {
        double[] probs = new double[distinctTrees.size()];
        Thread[] myThreads = new Thread[_numThread];

        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        gtp.setParallel(true);
        gtp.preProcess(speciesNetwork, distinctTrees, true);


        for(int i=0; i<_numThread; i++){
            myThreads[i] = new MyThreadFromScratch(gtp, speciesNetwork, distinctTrees, species2alleles, probs);
            myThreads[i].start();
        }

        for(int i=0; i<_numThread; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        /*
        for(int i=0; i<probs.length; i++){
            System.out.println(distinctTrees.get(i));
            System.out.println(probs[i]);
            System.out.println();
        }
        */
        double initialProb = 0;
        for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double totalProb = 0;
            for(int id: triple.Item3){
                totalProb += probs[id];
            }
            initialProb += Math.log(totalProb) * triple.Item2;
        }

        /*
        System.out.println("\n" + initialProb);
        System.out.println(System.currentTimeMillis()-start);
        System.exit(0);
        */

        return initialProb;
    }


    protected double computeProbability(Network speciesNetwork, List<Tree> geneTrees, Map<String, List<String>> species2alleles, List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList) {
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        double[] probList = new double[geneTrees.size()];
        gtp.calculateGTDistribution(speciesNetwork, geneTrees, species2alleles, probList);
        double total = 0;
        for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double subTotal = 0;
            for(int id: triple.Item3){
                subTotal += probList[id];
            }
            total += Math.log(subTotal) * triple.Item2;
        }
        return total;
    }




    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList, NetNode child, NetNode parent) {
        Set<NetNode> childNodes = new HashSet<NetNode>();
        childNodes.add(child);
        Set<NetNode> parentNodes = new HashSet<NetNode>();
        parentNodes.add(parent);
        if(_numThread>1){
            return computeProbabilityParallel(speciesNetwork, geneTrees, nbTreeAndCountAndBinaryIDList, childNodes, parentNodes);
        }
        else{
            return computeProbability(speciesNetwork, geneTrees, nbTreeAndCountAndBinaryIDList, childNodes, parentNodes);
        }    }


    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList, NetNode node, boolean changeBranchLength) {
        Set<NetNode> childNodes = new HashSet<NetNode>();
        Set<NetNode> parentNodes = new HashSet<NetNode>();
        if(changeBranchLength){
            //parentNodes = new HashSet<NetNode>();
            parentNodes.add(node);
            for(Object childObject: node.getChildren()){
                childNodes.add((NetNode)childObject);
            }

            childNodes.add(node);
            for(Object parentObject: node.getParents()){
                parentNodes.add((NetNode)parentObject);
            }
        }
        else{
            for(Object parentObject: node.getParents()){
                parentNodes.add((NetNode)parentObject);
            }
            childNodes.add(node);
        }

        if(_numThread>1){
            return computeProbabilityParallel(speciesNetwork, geneTrees, nbTreeAndCountAndBinaryIDList, childNodes, parentNodes);
        }
        else{
            return computeProbability(speciesNetwork, geneTrees, nbTreeAndCountAndBinaryIDList, childNodes, parentNodes);
        }
    }


    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList, Set<NetNode> childNodes, Set<NetNode> parentNodes) {
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        double[] probList = new double[geneTrees.size()];
        gtp.calculateGTDistribution(speciesNetwork, geneTrees, childNodes, parentNodes, probList);
        double total = 0;
        for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double subTotal = 0;
            for(int id: triple.Item3){
                subTotal += probList[id];
            }
            total += Math.log(subTotal) * triple.Item2;
        }
        return total;
    }


    public double computeProbabilityParallel(Network<Object> speciesNetwork, List<Tree> distinctTrees, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList, Set<NetNode> childNodes, Set<NetNode> parentNodes) {
        double[] probs = new double[distinctTrees.size()];
        Thread[] myThreads = new Thread[_numThread];

        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        gtp.setParallel(true);
        gtp.preProcess(speciesNetwork, distinctTrees, false);

        //System.out.println("\ngts:" +distinctTrees);
        for(int i=0; i<_numThread; i++){
            myThreads[i] = new MyThreadFromNonScratch(gtp, speciesNetwork, distinctTrees, childNodes, parentNodes, probs);
            myThreads[i].start();
        }

        for(int i=0; i<_numThread; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        double initialProb = 0;
        for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double totalProb = 0;
            for(int id: triple.Item3){
                totalProb += probs[id];
            }
            initialProb += Math.log(totalProb) * triple.Item2;
        }

        /*
        System.out.println("\n" + initialProb);
        System.out.println(System.currentTimeMillis()-start);
        System.exit(0);
        */

        return initialProb;
    }

    protected String network2String(Network speciesNetwork){
        RnNewickPrinter<Double> rnNewickPrinter = new RnNewickPrinter<Double>();
        StringWriter sw = new StringWriter();
        rnNewickPrinter.print(speciesNetwork, sw);
        return sw.toString();
    }

    protected String network2String(final DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork){
        Func1<String,String> _getNetworkNodeLabel  = new Func1<String,String>()
        {
            public String execute(String node) {
                return node;
            }
        };

        Func1<String, Iterable<String>> _getDestinationNodes = new Func1<String, Iterable<String>>() {
            public Iterable<String> execute(String node) {
                return new GetDirectSuccessors<String,PhyloEdge<String>>().execute(speciesNetwork, node);
            }
        };


        Func2<String, String, String> _getNetworkDistanceForPrint = new Func2<String, String, String>() {
            public String execute(String parent, String child) {
                PhyloEdge<String> edge = speciesNetwork.getEdge(parent, child);
                if(edge.getBranchLength()==null){
                    return null;
                }
                return edge.getBranchLength()+"";
            }
        };

        Func2<String, String, String> _getProbabilityForPrint = new Func2<String, String, String>() {
            public String execute(String parent, String child) {
                PhyloEdge<String> edge = speciesNetwork.getEdge(parent, child);
                if(edge.getProbability()==null){
                    return null;
                }
                return edge.getProbability()+"";
            }
        };

        Func2<String, String, String> _getSupportForPrint = new Func2<String, String, String>() {
            public String execute(String parent, String child) {
                PhyloEdge<String> edge = speciesNetwork.getEdge(parent, child);
                if(edge.getSupport()==null){
                    return null;
                }
                return edge.getSupport()+"";
            }
        };

        Func1<String, HybridNodeType> _getHybridTypeForPrint = new Func1<String, HybridNodeType>()
        {
            public HybridNodeType execute(String node)
            {
                int inDegree = new GetInDegree<String,PhyloEdge<String>>().execute(speciesNetwork, node);
                return inDegree == 2 ? HybridNodeType.Hybridization : null;
            }
        };

        Func1<String,String> _getHybridNodeIndexForPrint = new Func1<String, String>() {
            List<String> hybridNodes = new ArrayList<String>();

            public String execute(String node) {
                int inDegree = new GetInDegree<String,PhyloEdge<String>>().execute(speciesNetwork, node);
                if(inDegree == 2){
                    int index = hybridNodes.indexOf(node) + 1;
                    if(index == 0){
                        hybridNodes.add(node);
                        return hybridNodes.size()+"";
                    }
                    else{
                        return index + "";
                    }
                }
                else{
                    return null;
                }
            }
        };

        try{
            StringWriter sw = new StringWriter();
            //   new RichNewickPrinterCompact<String>().print(true, "R", _getNetworkNodeLabel, _getDestinationNodes, _getNetworkDistanceForPrint, _getSupportForPrint, _getProbabilityForPrint, _getHybridNodeIndexForPrint, _getHybridTypeForPrint, sw);
            RichNewickPrinterCompact<String> printer = new RichNewickPrinterCompact<String>();
            printer.setGetBranchLength(_getNetworkDistanceForPrint);
            printer.setGetProbability(_getProbabilityForPrint);
            printer.setGetSupport(_getSupportForPrint);

            printer.print(true, new FindRoot<String>().execute(speciesNetwork), _getNetworkNodeLabel, _getDestinationNodes, _getHybridNodeIndexForPrint, _getHybridTypeForPrint, sw);
            sw.flush();
            sw.close();
            return sw.toString();
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return null;
    }




    protected Network networkNew2Old(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork){
        Network<Object> bniNetwork = null;
        try{

            String networkString = network2String(speciesNetwork);
            bniNetwork = string2Network(networkString);
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return bniNetwork;
    }

    protected Network string2Network(String networkString){
        try{
            RichNewickReaderAST reader = new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
            reader.setHybridSumTolerance(BigDecimal.valueOf(0.00001));
            NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            RichNewickReadResult<Networks> readResult = reader.read(new ByteArrayInputStream(networkString.getBytes()));
            return transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return null;
    }

}
