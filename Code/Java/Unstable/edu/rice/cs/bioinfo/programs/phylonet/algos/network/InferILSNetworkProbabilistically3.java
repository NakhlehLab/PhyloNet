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
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkWholeNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours.AllNeighboursHillClimberSteepestAscent;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
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
import java.lang.Math;

/**
 * Created with IntelliJ IDEA.
 * User: jianrongdong
 * Date: 3/18/13
 * Time: 9:45 AM
 * To change this template use File | Settings | File Templates.
 */

/* This is for cross validation. It is modified from Yun's InferILSNetworkProbabilistically */

public class InferILSNetworkProbabilistically3 extends MDCOnNetworkYFFromRichNewickJung {
    private Network[] _optimalNetworks;
    private double[] _optimalScores;
    int _maxRounds;
    int _maxTryPerBranch;
    double _improvementThreshold;
    double _maxBranchLength;
    double _Brent1;
    double _Brent2;
    int _numFolds = 10; // default as 10-folds stuff.

    // No change needed
    public InferILSNetworkProbabilistically3(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }

    // Another constructor added to input number of folds in cross validation
    public InferILSNetworkProbabilistically3(int numFolds) {
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
        this._numFolds = numFolds;
    }

    // No change needed
    public void setBrentParameter(int maxRounds,
                                  int maxTryPerBranch,
                                  double improvementThreshold,
                                  double maxBranchLength,
                                  double Brent1,
                                  double Brent2){
        _maxRounds = maxRounds;
        _maxTryPerBranch = maxTryPerBranch;
        _improvementThreshold = improvementThreshold;
        _maxBranchLength = maxBranchLength;
        _Brent1 = Brent1;
        _Brent2 = Brent2;
    }

    // CV, changed from inferNetwork
    public long CV(List<Tree> gts,Map<String,List<String>> species2alleles,Long maxExaminations,
                   Long maxReticulations, int diameterLimit, Network startNetwork, int numSol){
        _optimalNetworks = new Network[numSol];
        _optimalScores = new double[numSol];

        List<Tree>gts2 = new ArrayList<Tree>(gts);

        List<Tree> distinctTrees = new ArrayList<Tree>();

        List<Tuple3<Tree, List<Double>, List<Integer>>> nbTreeAndCountAndBinaryIDListForCV =
                new ArrayList<Tuple3<Tree, List<Double>, List<Integer>>>();  // New data structure

        summarizeGeneTreestoKfolds(gts2, distinctTrees, nbTreeAndCountAndBinaryIDListForCV);   // Replace old routine

        /*  1. The outside loop is for the number of reticulation edges needed. Since it is undetermined, I need to
               go up until data analysis shows to stop, or it exceeds the maximum reticulation edges limit.

            2. The inside loop is for cross validation. K rounds. No model with all data.

            3. MyContainer holds the tsublist and vsublist MSE results for every sublist and every
               reticulation node. MyContainer2 holds the total training MSE and total validation
               MSE results for every reticulation node.

            4: The internal process inside the double loop.

            a. Given certain k (reticulation nodes), build a model with the tsublist.
               Return ML value, probList (for each distinct gene tree in the whole list),
               and the optimized network.
            b. Save the maximum likelihood into MyContainer.
            c. Use the probList to compute the tsublist MSE.
            d. Save the tsublist MSE squared into MyContainer.
            e. Use the optimized network and the sub-list to compute the validation total errors
            f. Save the validation sub-list MSE squared into MyContainer.

            5. If the weight of a tree is zero in a sublist, do not pass tuple information to the
            sublist.

            6. The external process outside the double loop.
               MSE_vsublist = \sum{vsublist error squared}/numDistinctTrees in vsublist.
               MSE_tsublist = \sum{tsublist error squared}/numDistinctTrees in tsublist;
         */

        // External loop to increase the reticulation nodes.
        boolean isTurningPointFound = false;
        long curMaxReticulations = 0;  // initial value of reticulation nodes.
        long correctK = -1; // the correct number of reticulation nodes

        ArrayList<Double> myContainer_t = new ArrayList<Double>();  // container for tsublists (1 to K) MSE for k=i
        ArrayList<Double> myContainer_v = new ArrayList<Double>();  // container for vsublists (1 to K) MSE for k=i
        ArrayList<Double> myContainer2_t = new ArrayList<Double>(); // container for total MSE for all reticulations
        ArrayList<Double> myContainer2_v = new ArrayList<Double>(); // container for total MSE for all reticulations

        while (!isTurningPointFound){
            System.out.println(); // restart from a line

            // Loop for K tsublists to train the network models
            for (int s=_numFolds+1; s<= _numFolds*2; s++) {
                System.out.println("curMaxReticulations = " + curMaxReticulations+ "  s = "+s);

                // Transfer info from nbTreeAndCountAndBinaryIDListForCV to a new
                // nbTreeAndCountAndBinaryIDList used for tsublist
                List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList =
                        new ArrayList<Tuple3<Tree, Double, List<Integer>>>();
                for (Tuple3<Tree,List<Double>,List<Integer>> tuple: nbTreeAndCountAndBinaryIDListForCV) {
                    if (tuple.Item2.get(s) > 0) { // weight > 0 --- zero weight terms not included
                        Tuple3<Tree,Double,List<Integer>> triple =
                                new Tuple3<Tree,Double,List<Integer>>(tuple.Item1,
                                                                      tuple.Item2.get(s),
                                                                      tuple.Item3);
                        nbTreeAndCountAndBinaryIDList.add(triple);
                    }
                }

                Arrays.fill(_optimalScores, Double.NEGATIVE_INFINITY);

                // Use the correct part of the gts to build the model.
                int numVsublistTrees = gts2.size()/_numFolds;
                List<Tree>gts3 = new ArrayList<Tree>();
                int s1 = s-_numFolds; // the section of gts not used --- _numFolds*(s1-1) to _numFolds*s1 (exclusive)
                if (s1 == 1) {
                    gts3.addAll(gts2.subList(numVsublistTrees,gts2.size()));
                }
                else {
                    gts3.addAll(gts2.subList(0,numVsublistTrees*(s1-1)));
                    gts3.addAll(gts2.subList(numVsublistTrees*s1,gts2.size()));
                }

                DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork =
                        getStartNetwork(gts3, species2alleles, startNetwork);  // changed from gts to gts3

                NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,
                        String,PhyloEdge<String>> allNeighboursStrategy =
                        new NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,
                                String, PhyloEdge<String>>(makeNode, makeEdge);

                AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,
                        String,PhyloEdge<String>,Double> searcher =
                        new AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,
                                String, PhyloEdge<String>, Double>(allNeighboursStrategy);

                Comparator<Double> comparator = getDoubleScoreComparator();

                Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> scorer =
                        getScoreFunction(distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);

                HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result =
                        searcher.search(speciesNetwork, scorer, comparator, maxExaminations,
                                        curMaxReticulations, diameterLimit);
                // Results are saved in _optimalNetworks and _optimalScores. Choose the best ones.
                Network curSpeciesNetwork = _optimalNetworks[0];
                Double curBestScore = _optimalScores[0];


                // Get theoretical frequency of each distinct gene tree with tsublist trained network
                GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
                List<Double> probList;
                // The computed probability list is for the distinct binary trees in the whole list.
                probList = gtp.calculateGTDistribution(curSpeciesNetwork,
                                                       distinctTrees,
                                                       species2alleles);

                // Compute real frequency for each distinct binary tree (distinctTrees)
                // in tsublist and vsublist respectively.
                //
                // The total weight of a non-binary tree, its binary refinements,
                // and the theoretical probability of each distinct binary tree are known.
                // Pick the binary trees with highest theoretical probability and put
                // the averaged weight on them.

                // 1. Work on nbTreeAndCountAndBinaryIDList (tsublist).
                List<Double> realWeightsDistinctTrees_tsublist = new ArrayList<Double>();
                for (int i=0; i<=distinctTrees.size()-1;i++)
                    realWeightsDistinctTrees_tsublist.add(0.0); // initialize weights to 0.0

                for (Tuple3<Tree,Double,List<Integer>> triple: nbTreeAndCountAndBinaryIDList) {
                    // For each non-binary tree, find the set of binary gene trees
                    // with the highest probability values.

                    // Get maxProb value
                    double maxProb = 0.0;
                    for(int id: triple.Item3){ // id is the index in discreteTrees and in probList
                        maxProb = Math.max(maxProb, probList.get(id));  // get maxProb value
                    }

                    // Get the number of corresponding binary trees in distinctTrees that have
                    // maxProb and their id's in discreteTrees.
                    int maxProbCount = 0;
                    List<Integer> maxProbBinaryIDList = new ArrayList<Integer>();
                    for (int id: triple.Item3) {
                        if (Math.abs(maxProb - probList.get(id))/maxProb<1e-10) {
                            // Essentially maxProb = probList.get(id)
                            maxProbCount++;            // number of trees have maxProb
                            maxProbBinaryIDList.add(id);   // the tree ids that have maxProb
                        }
                    }

                    // Compute the average weight.
                    double weight = triple.Item2;  // get the non-binary tree weight
                    double averageWeight = weight / (double) maxProbCount;  // average weight to the GTs.

                    // Add weight to realWeightsDistrinctTrees_tsublist
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
                ArrayList<Double> realFreqsDistinctTrees_tsublist =
                        new ArrayList<Double>();
                for (int i=0; i<=realWeightsDistinctTrees_tsublist.size()-1;i++)
                    realFreqsDistinctTrees_tsublist.add(realWeightsDistinctTrees_tsublist.get(i)/sum);
                // Compute the MSE for the tsublist
                double SE = 0;
                double oneError;
                for (int i=0; i<=realFreqsDistinctTrees_tsublist.size()-1;i++) {
                    if (realFreqsDistinctTrees_tsublist.get(i)>0.0) { // consider only types existent in tsublist
                        oneError = realFreqsDistinctTrees_tsublist.get(i)-probList.get(i);
                        SE += oneError*oneError;
                    }
                }
                double MSE = SE/count;

                // Store MSE into myContainer
                myContainer_t.add(MSE);


                // 2. Work on nbTreeAndCountAndBinaryIDList_vsublist.
                // Transfer info from nbTreeAndCountAndBinaryIDListForCV
                List<Tuple3<Tree,Double, List<Integer>>> nbTreeAndCountAndBinaryIDList_vsubList =
                        new ArrayList<Tuple3<Tree,Double, List<Integer>>>();

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
                for (int i = 0; i <= distinctTrees.size()-1;i++)
                    realWeightsDistinctTrees_vsublist.add(0.0); // initialize weights to 0.0


                for (Tuple3<Tree,Double,List<Integer>> triple: nbTreeAndCountAndBinaryIDList_vsubList) {
                    // For each non-binary tree, find the set of binary gene trees
                    // with the highest probability values.

                    // Get maxProb value
                    double maxProb = 0.0;
                    for (int id: triple.Item3){ // id is the index in discreteTrees and in probList
                        maxProb = Math.max(maxProb, probList.get(id)); // get maxProb value
                    }

                    // Get the number of corresponding binary trees in distinctTrees that have
                    // maxProb and their id's in discreteTrees.
                    int maxProbCount = 0;
                    List<Integer> maxProbBinaryIDList = new ArrayList<Integer>();
                    for (int id: triple.Item3) {
                        if (Math.abs(maxProb - probList.get(id))/maxProb<1e-10) {
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

                // Compute the MSE for the vsublist
                SE = 0;
                for (int i=0; i<=realFreqsDistinctTrees_vsublist.size()-1;i++)
                    if (realFreqsDistinctTrees_vsublist.get(i)>0.0) {
                        oneError = realFreqsDistinctTrees_vsublist.get(i)-probList.get(i);
                        SE += oneError*oneError;
                    }
                MSE = SE/count;

                // Store MSE into myContainer
                myContainer_v.add(MSE);

                // Note: Since we are only using a sub-list, it may not have
                // all the distinct gene trees as in the distinctGTs. Therefore,
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

            for (double entry: myContainer2_t) {
                System.out.println("In myContainer2_t:" + entry);
            }
            for (double entry: myContainer2_v) {
                System.out.println("In myContainer2_v:" + entry);
            }

            if (myContainer2_t.size() <= 1) { // get two data points before any analysis
                curMaxReticulations++;
                continue;   // skip the rest of the code and go for the while loop
            }
            else { // >= 2 data points

                // If MSE_t not always going down, something wrong. Stop.
                for (int i=0; i<=myContainer2_t.size()-2;i++) {
                    double MSE1_t = myContainer2_t.get(i);
                    double MSE2_t = myContainer2_t.get(i+1);
                    if (MSE2_t >= MSE1_t) {
                        System.out.println("Unbelievable! MSE even increase in the training set.");
                        System.exit(-1);
                    }
                }

                // There are only two data points and goes up at the very beginning -- tree is the answer.
                double MSE0_v = myContainer2_v.get(0);
                double MSE1_v = myContainer2_v.get(1);
                if (myContainer2_v.size() == 2 && MSE0_v < MSE1_v) {
                    System.out.println("Tree is the best solution. Number of reticulation nodes = 0");
                    System.exit(0);
                }

                // Need more data
                if (myContainer2_v.size()==2) {  // Only 0 and 1 reticulation nodes but tree is not best, is not enough.
                    curMaxReticulations++;
                    continue; // skip the rest of the code in the while loop
                }

                // Now there are >= 3 MSE results.
                // Loop
                //      if there is a turning point,
                //          Set the best reticulation node number.
                //          Jump out of the while loop (break) and skip the next if (set isEnough to true)

                for (int i=0; i <= myContainer2_v.size()-3;i++) {

                    double MSE1 = myContainer2_v.get(i);
                    double MSE2 = myContainer2_v.get(i+1);
                    double MSE3 = myContainer2_v.get(i+2);
                    if (MSE2 < MSE1 && MSE2 < MSE3) { // turn point found
                        System.out.println("Turning point found. The right reticulation number is " + (i+1));
                        correctK = (long)(i+1);
                        isTurningPointFound = true;
                        break;
                    }
                }

                // No turning point
                if (!isTurningPointFound)
                    if (myContainer2_v.size() == maxReticulations+2) { // Already tested 0 to maxReticulation+1 nodes
                        System.out.println("Tested 0, 1, ..., maxReticulations, and maxReticulation + 1 nodes. No turning point found.");
                        System.exit(-1);
                    }
                    else
                        curMaxReticulations++;
            }  // else >=2 data points
        } // while
        for (double entry: myContainer2_t) {
            System.out.println("In myContainer2_t:" + entry);
        }
        for (double entry: myContainer2_v) {
            System.out.println("In myContainer2_v:" + entry);
        }
        return correctK;
    }

    // No change needed
    private DirectedGraphToGraphAdapter<String,PhyloEdge<String>> getStartNetwork(List<Tree> gts, Map<String,List<String>> species2alleles, Network<Object> startingNetwork){
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
                sol = mdc.inferSpeciesTree(gts, false, 1, false, 100, true, -1).get(0);
            }
            else{
                sol = mdc.inferSpeciesTree(gts, allele2species, false, 1, false, 100, true, -1).get(0);
            }

            Tree startingTree= Trees.generateRandomBinaryResolution(sol._st);
            startingNetwork = string2Network(startingTree.toString());
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
        }

        String newNetwork = network2String(startingNetwork);

        return makeNetwork(newNetwork);

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

        // Check types of gene trees and find the number of each type in originalGTs
        for(Tree tr: originalGTs){

            // The structure of nbTreeAndCountAndBinaryIDList. A triple = (nbTree, CountList, BinaryIDList)
            //
            // nbTree is a gene tree, it could be non-binary.  In original program,
            //
            // Count = the number of effective occurrences (all weights added together) of the nbTree in
            // originalGTs.
            //
            // Expand Count to CountList in nbTreeAndCountAndBinaryIDListForCV into an array of 2K+1 integers.
            //
            // The 1st double = effective occurrence count in the originalOriginalGTs.
            // The 2nd double = effective occurrence count in the 1stSubListOriginalGTs.
            // ...
            // The (K+1)th double = effective occurrence count in the (K+1)thSubListOriginalGTs.
            // The (K+1+1)th double = effective occurrence count in the 1stComplementarySubListOriginalGTs
            // The (K+1+2)th double = effective occurrence count in the 2ndComplementarySubListOriginalGTs
            // ...
            // The (2K+1)th double = effective occurrence count in the KthComplementarySubListOriginalGTs
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

            double weight = 1.0;
            // double weight = ((STINode<Double>)tr.getRoot()).getData();
            int index = 0;     // the index of nbTreeAndCountAndBinaryIDList
            Tuple3<Tree, List<Double>, List<Integer>> newTuple = null;  // initial tuple

            for(Tuple3<Tree, List<Double>, List<Integer>> triple: nbTreeAndCountAndBinaryIDListForCV){
                // loop for all tuples in nbTreeAndCountAndBinaryIDList up to now.

                if(Trees.haveSameRootedTopology(tr, triple.Item1)){
                    // The same topology tree, tr matches a triple's tree
                    // in nbTreeAndCountAndBinaryIDList.
                    //
                    // Update the triple entry and give it to newTuple, and jump out of this for loop
                    //
                    triple.Item2.set(0, triple.Item2.get(0)+weight);  // update the total occurrence stored in Item2
                    newTuple = new Tuple3<Tree, List<Double>, List<Integer>>(tr, triple.Item2, triple.Item3);
                    // Create a newTuple = an updated version of the one saved in
                    // nbTreeAndCountAndBinaryIDList

                    break;  // when breaking out, the index points to this item
                }
                index++;  // if no match, increment index to be the next item.

            } // for loop for all tuples in nbTreeAndCountAndBinaryIDList up to now.

            if(newTuple!=null){ // tr Matched triple and jumped out case,
                // put updated newTuple into nbTreeAndCountAndBinaryIDList at index position
                // that is, to update nbTreeAndCountAndBinaryIDList

                nbTreeAndCountAndBinaryIDListForCV.set(index, newTuple);
            }
            else{  // Not match case, tr did not match a triple's Item 1 in nbTreeAndCountAndBinaryIDList.
                // Create a new binaryID for each of tr's binary resolutions and
                // put them into binaryIDs.

                List<Integer> binaryIDs = new ArrayList<Integer>();
                for(Tree btr: Trees.getAllBinaryResolution(tr)){
                    // loop for all binary resolution of tr, tr could be non-binary
                    // btr is one of its binary resolution trees.
                    // distinctGTs are all binary trees.

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

                // Create a new Tuple to be put into nbTreeAndCountAndBinaryIDList
                // Create the Countlist
                ArrayList<Double> temp = new ArrayList<Double>(_numFolds*2+1);
                temp.add(0,weight);  // the first element of the CountList
                for (int i=1; i<=_numFolds*2;i++)  // the other elements of the CountList
                    temp.add(0.0);

                newTuple = new Tuple3<Tree, List<Double>, List<Integer>>(tr, temp, binaryIDs);

                nbTreeAndCountAndBinaryIDListForCV.add(newTuple);
                // This new tr is a first encounter. Together with its binaryID's in
                // distinctGT, a newTuple is created and is put into
                // nbTreeAndCountAndBinaryIDList
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

        // Loop for all vsublists, compute the correct frequencies
        for (int i=0; i<=_numFolds-1;i++) {
            List<Tree> aSubListOriginalGTs = originalGTs.subList(numASubListOriginalGTs*i,
                                                                 numASubListOriginalGTs*(i+1));  // create vsublist

            for(Tree tr: aSubListOriginalGTs){ // loop for the sub-list of original gene trees.
                // double weight = ((STINode<Double>)tr.getRoot()).getData();
                double weight = 1.0;
                int index = 0;     // the index of nbTreeAndCountAndBinaryIDList
                Tuple3<Tree, List<Double>, List<Integer>> newTuple = null;  // initial tuple

                for(Tuple3<Tree, List<Double>, List<Integer>> triple: nbTreeAndCountAndBinaryIDListForCV){
                    // Loop for all tuples in nbTreeAndCountAndBinaryIDList.

                    if(Trees.haveSameRootedTopology(tr, triple.Item1)){
                        // The same topology tree, tr matches a triple's tree in nbTreeAndCountAndBinaryIDList.
                        // Update CountList and give it to newTuple, and jump out of this loop

                        double temp1 = triple.Item2.get(i+1);  // Get the original value of the i+1 th term
                                                            // in CountList
                        temp1 +=weight;  // increment temp1
                        triple.Item2.set(i+1,temp1); // increment the (i+1)th element in CountList

                        newTuple = new Tuple3<Tree, List<Double>, List<Integer>>(tr, triple.Item2, triple.Item3);
                        // create a newTuple = an updated version of the one saved in nbTreeAndCountAndBinaryIDList

                        break;  // when breaking out, the index points to this item
                    }
                    index++;  // if no match, increment index to be the next item.

                } // for loop for all tuples in nbTreeAndCountAndBinaryIDList.
                // Note: In this round every tree must find its counterpart in
                // nbTreeAndCountAndBinaryIDList. Otherwise throw an exception.

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

    // No changes needed
    private Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }

    // No changes needed
    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                //System.out.println("Start scoring ...");
                //System.out.println(network2String(network));
                //long start = System.currentTimeMillis();
                Network<Object> speciesNetwork = networkNew2Old(network);

                double score = findOptimalBranchLength(speciesNetwork, distinctTrees, species2alleles,nbTreeAndCountAndBinaryIDList);

                if(score > _optimalScores[_optimalNetworks.length-1]){
                    boolean exist = false;
                    for(int i=0; i<_optimalNetworks.length; i++){
                        if(_optimalNetworks[i]==null)break;

                        if(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeClusterDistance(speciesNetwork, _optimalNetworks[i])[2]<0.000001 &&
                                edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeTripartitionDistance(speciesNetwork, _optimalNetworks[i])[2]<0.000001){
                            exist = true;
                            break;
                        }
                    }
                    if(!exist){
                        int index = -1;
                        for(int i=0; i<_optimalScores.length; i++){
                            if(score > _optimalScores[i]){
                                index = i;
                                break;
                            }
                        }
                        for(int i=_optimalScores.length-1; i>index; i--){
                            _optimalNetworks[i] = _optimalNetworks[i-1];
                            _optimalScores[i] = _optimalScores[i-1];
                        }
                        _optimalScores[index] = score;
                        _optimalNetworks[index] = string2Network(network2String(speciesNetwork));
                        //System.out.println(network2String(speciesNetwork) + ": "+score);
                    }
                }
                //System.out.println("End scoring ..." + (System.currentTimeMillis()-start)/1000.0);
                return score;
            }
        };
    }

    // No changes needed
    private double findOptimalBranchLength(final Network<Object> speciesNetwork, final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList){
        boolean continueRounds = true; // keep trying to improve network

        //boolean isNetwork = false;
        for(NetNode<Object> node: speciesNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent,1.0);
                if(node.isNetworkNode()){
                    node.setParentProbability(parent, 0.5);
                    //isNetwork = true;
                }
            }
        }

        //final boolean hasNetworkNode = isNetwork;
        //long start = System.currentTimeMillis();
        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList));  // records the GTProb of the network at all times
        //System.out.println("first one : " + (System.currentTimeMillis()-start)/1000.0);
        //final Container<Double> timeTotal = new Container<Double>(System.currentTimeMillis()/1000.0);  // records the GTProb of the network at all times
        //final Container<Integer> timeCount = new Container<Integer>(0);  // records the GTProb of the network at all times


        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
            /*
            * Prepare a random ordering of network edge examinations each of which attempts to change a branch length or hybrid prob to improve the GTProb score.
            */

            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();
            List<Proc> assigmentActions = new ArrayList<Proc>(); // store adjustment commands here.  Will execute them one by one later.

            // add branch length adjustments to the list

            for(final NetNode<Object> parent : speciesNetwork.bfs())
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
                            //timeTotal.setContents(System.currentTimeMillis()/1000.0);
                            //final Container<Integer> callCount = new Container<Integer>(0);
                            UnivariateFunction functionToOptimize = new UnivariateFunction() {
                                public double value(double suggestedBranchLength) {  // brent suggests a new branch length
                                    //callCount.setContents(callCount.getContents()+1);
                                    double incumbentBranchLength = child.getParentDistance(parent);


                                    // mutate and see if it yields an improved network


                                    child.setParentDistance(parent, suggestedBranchLength);

                                    double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, parent);
/*
                                            timeCount.setContents(timeCount.getContents()+1);
                                            if(timeCount.getContents()%1000==0){
                                                //timeTotal.setContents(timeTotal.getContents()+(System.currentTimeMillis()/1000.0-timeTotal.getContents()));
                                                System.out.println(timeCount.getContents() + ": " +(System.currentTimeMillis()/1000.0-timeTotal.getContents())/timeCount.getContents());
                                            }
*/
                                    /*
                                    RnNewickPrinter<Object> rnNewickPrinter = new RnNewickPrinter<Object>();
                                    StringWriter sw = new StringWriter();
                                    rnNewickPrinter.print(speciesNetwork, sw);
                                       String inferredNetwork = sw.toString();
                                        System.out.println(inferredNetwork + "\t" + lnProb);
                                    */
                                    if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                                    {
                                        lnGtProbOfSpeciesNetwork.setContents(lnProb); // System.out.println("(improved)");
                                        //System.out.println(network2String(speciesNetwork) + ":" + lnProb);
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

                            //System.out.println(computeProbability(speciesNetwork, geneTrees, counter, child, parent) + " vs. " + computeProbability(speciesNetwork, geneTrees, counter));

                            //   System.out.println("-----------------------------------------------------------------------");
                            //System.out.println(callCount.getContents());
                            computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, parent);
                        }
                    });
                }
            }


            // add hybrid probs to hybrid edges
            for(final NetNode<Object> child : speciesNetwork.getNetworkNodes()) // find every hybrid node
            {

                Iterator<NetNode<Object>> hybridParents = child.getParents().iterator();
                final NetNode hybridParent1 = hybridParents.next();
                final NetNode hybridParent2 = hybridParents.next();
                //child.setParentProbability(hybridParent1,0.5);
                //child.setParentProbability(hybridParent2,0.5);
                assigmentActions.add(new Proc()
                {
                    public void execute()
                    {
                        UnivariateFunction functionToOptimize = new UnivariateFunction() {
                            public double value(double suggestedProb) {

                                double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);

                                // try new pair of hybrid probs
                                child.setParentProbability(hybridParent1, suggestedProb);
                                child.setParentProbability(hybridParent2, 1.0 - suggestedProb);

                                //long start = System.currentTimeMillis();
                                double lnProb = computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, null);
                                //System.out.print((System.currentTimeMillis()-start)/1000.0 + " ");
                                //double lnProb = computeProbability(speciesNetwork, geneTrees, counter);
                                /*
                                timeCount.setContents(timeCount.getContents()+1);
                                if(timeCount.getContents()%1000==0){
                                    //timeTotal.setContents(timeTotal.getContents()+(System.currentTimeMillis()/1000.0-timeTotal.getContents()));
                                    System.out.println(timeCount.getContents() + ": " +(System.currentTimeMillis()/1000.0-timeTotal.getContents())/timeCount.getContents());
                                }
                                */
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
                        computeProbability(speciesNetwork, distinctTrees, nbTreeAndCountAndBinaryIDList, child, null);
                        //System.out.println(computeProbability(speciesNetwork, geneTrees, counter, child, null) + " vs. " + computeProbability(speciesNetwork, geneTrees, counter));

                    }
                });


            }

            //     Collections.shuffle(assigmentActions); // randomize the order we will try to adjust network edge properties

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

        //System.out.println(standardComputeProbability(speciesNetwork, geneTrees, counter) + " vs. " + lnGtProbOfSpeciesNetwork.getContents());
        return lnGtProbOfSpeciesNetwork.getContents();
    }





    // No changes
    private double computeProbability(Network speciesNetwork, List<Tree> geneTrees, Map<String, List<String>> species2alleles, List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList) {
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        List<Double> probList = gtp.calculateGTDistribution(speciesNetwork, geneTrees, species2alleles);
        double total = 0;
        for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double maxProb = 0;
            for(int id: triple.Item3){
                maxProb = Math.max(maxProb, probList.get(id));
            }
            total += Math.log(maxProb) * triple.Item2;
        }
        return total;
    }

    // No changes
    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple3<Tree, Double, List<Integer>>> nbTreeAndCountAndBinaryIDList, NetNode child, NetNode parent) {
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        List<Double> probList = gtp.calculateGTDistribution(speciesNetwork, geneTrees, child, parent);
        double total = 0;
        for(Tuple3<Tree, Double, List<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double maxProb = 0;
            for(int id: triple.Item3){
                maxProb = Math.max(maxProb, probList.get(id));
            }
            total += Math.log(maxProb) * triple.Item2;
        }
        return total;
    }

    // No changes
    private String network2String(Network speciesNetwork){
        RnNewickPrinter<Double> rnNewickPrinter = new RnNewickPrinter<Double>();
        StringWriter sw = new StringWriter();
        rnNewickPrinter.print(speciesNetwork, sw);
        return sw.toString();
    }

    // No changes
    private String network2String(final DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork){
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



    // No changes
    private Network networkNew2Old(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork){
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

    // No changes
    private Network string2Network(String networkString){
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