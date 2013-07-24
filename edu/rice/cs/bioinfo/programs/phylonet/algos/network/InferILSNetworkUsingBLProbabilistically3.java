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
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BfsSearch;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeClusterWD;
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
 * Date: 5/1/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferILSNetworkUsingBLProbabilistically3 extends MDCOnNetworkYFFromRichNewickJung {
    private Network[] _optimalNetworks;
    private double[] _optimalScores;
    private int _maxRounds;
    private int _maxTryPerBranch;
    private double _improvementThreshold;
    private double _maxBranchLength;
    private double _Brent1;
    private double _Brent2;
    private Long _maxExaminations;
    private int _diameterLimit;
    private Network _startNetwork;
    protected int _numFolds = 10; // default as 10-folds cross-validation. Can be changed with non-default constructor.

    private Map<SpeciesPair, Double> _pair2time;


    public InferILSNetworkUsingBLProbabilistically3(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }

    public void setSearchParameter(int maxRounds,
                                   int maxTryPerBranch,
                                   double improvementThreshold,
                                   double maxBranchLength,
                                   double Brent1,
                                   double Brent2,
                                   Long maxExaminations,
                                   int diameterLimit,
                                   Network startNetwork){
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
    }

    public int CV(List<Tree> gts, Map<String,List<String>> species2alleles, int maxReticulations, int numSol){
        _optimalNetworks = new Network[numSol];
        _optimalScores = new double[numSol];
        Arrays.fill(_optimalScores, Double.NEGATIVE_INFINITY);

        Set<String> speciesSet = new LinkedHashSet<String>();
        Map<String, String> allele2species = computeAlleleToSpecies(gts, species2alleles, speciesSet);

        String[] snTaxa = new String[speciesSet.size()];
        int index = 0;
        for(String species: speciesSet){
            snTaxa[index++] = species;
        }


        // Note: when branch length is used, treat every tree as a distinct tree.
        /*  This is the explanation of the following codes.

            1. The external loop is for the number of reticulations.  Undetermined, go up until stoppage criteria are
               met, or exceeding the maximum reticulations limit.

            2. The internal loop is for cross validation. K rounds for K-fold cross validation.

            3. myContainer_t and myContainer_v hold the tsublist (training sublist) and vsublist (validation sublist)
               MSE results for every sublist and every reticulation node. myContainer2_t and myContainer2_v hold
               the total training MSE and total validation MSE results for every reticulation node.

            4: The process of the internal loop.
            a. Given certain k (reticulation nodes), build a model with the tsublist. Return ML value, probList,
               and the optimized network.
            b. Use probList to compute training MSE and save it into myContainer_t.
            c. Use probList to compute validation MSE and save it into myContainer_r.

            5. The external process outside the double loop.
               MSE_vsublist = \sum{vsublist error squared}/numTrees in vsublist.
               MSE_tsublist = \sum{tsublist error squared}/numTrees in tsublist;
        */

        // External loop to increase the reticulation nodes.
        boolean isTurningPointFound = false;
        int curMaxReticulations = 0;  // initial value of reticulation nodes.
        int correctK = -1; // the correct number of reticulation nodes

        List<Double> myContainer_t = new ArrayList<Double>();  // container for tsublists (1 to K) MSE for k=i
        List<Double> myContainer_v = new ArrayList<Double>();  // container for vsublists (1 to K) MSE for k=i
        List<Double> myContainer2_t = new ArrayList<Double>(); // container for total MSE for all reticulations
        List<Double> myContainer2_v = new ArrayList<Double>(); // container for total MSE for all reticulations

        // Add the likelihood counterparts
        List<Double> myContainer_t_ML = new ArrayList<Double>();  // container for tsublists (1 to K) ML for k=i
        List<Double> myContainer_v_ML = new ArrayList<Double>();  // container for vsublists (1 to K) ML for k=i
        List<Double> myContainer2_t_ML = new ArrayList<Double>(); // container for total ML for all reticulations
        List<Double> myContainer2_v_ML = new ArrayList<Double>(); // container for total ML for all reticulations


        while (!isTurningPointFound){
            System.out.println(); // restart from a line

            // Loop for K tsublists to train the network models
            for (int s=1; s<= _numFolds; s++) {
                System.out.println("curMaxReticulations = " + curMaxReticulations + "  s = "+s);

                Arrays.fill(_optimalScores, Double.NEGATIVE_INFINITY);

                // Put the part of gts for validation into gts2, and the part of gts to build the model into gts3.
                int numVsublistTrees = gts.size()/_numFolds;
                List<Tree> gts2 = new ArrayList<Tree>();
                List<Tree> gts3 = new ArrayList<Tree>();
                int s1 = s; // the section of gts not used --- _numFolds*(s1-1) to _numFolds*s1 (exclusive)
                if (s1 == 1) {
                    gts2.addAll(gts.subList(0,numVsublistTrees));
                    gts3.addAll(gts.subList(numVsublistTrees,gts.size()));
                }
                else {
                    gts2.addAll(gts.subList(numVsublistTrees*(s1-1),numVsublistTrees*s1));
                    gts3.addAll(gts.subList(0,numVsublistTrees*(s1-1)));
                    gts3.addAll(gts.subList(numVsublistTrees*s1,gts.size()));
                }


                computePairwiseCoalesceTime(gts3, allele2species, snTaxa);

                DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork =
                        getStartNetwork(gts3, allele2species,_startNetwork);

                NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,
                                                   String,
                                                   PhyloEdge<String>> allNeighboursStrategy =
                    new NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,
                                                           String,
                                                           PhyloEdge<String>>(makeNode, makeEdge);

                AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,
                                                       String,
                                                       PhyloEdge<String>,Double> searcher =
                    new AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,
                                                               String,
                                                               PhyloEdge<String>,
                                                               Double>(allNeighboursStrategy);

                //double[] operationProb = {0.15,0.15,0.2,0.5};
                //NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(operationProb, makeNode, makeEdge);
                //AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);

                Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> scorer =
                        getScoreFunction(gts3, species2alleles);

                Comparator<Double> comparator = getDoubleScoreComparator();

                //DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = getStartNetwork(gts, species2alleles,startNetwork);

                HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result =
                        searcher.search(speciesNetwork, scorer, comparator,
                                        _maxExaminations, maxReticulations, _diameterLimit); // search starts here

                //HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, maxExaminations, maxReticulations, new Long(100)); // search starts here

                // Results are saved in _optimalNetworks and _optimalScores. Choose the best ones.
                Network curSpeciesNetwork = _optimalNetworks[0];
                Double curBestScore = _optimalScores[0];

                // The following computes theoretical probabilities and real frequencies with the trained
                // network and the validation data set.

                // Get theoretical probability/frequency of each gene tree from the tsublist trained network
                GeneTreeWithBranchLengthProbabilityYF gtp = new GeneTreeWithBranchLengthProbabilityYF();
                List<Double> probList;
                // The theoretical probability list is for each tree.
                probList = gtp.calculateGTDistribution(curSpeciesNetwork, gts3, species2alleles);

                // Compute real frequency for each tree in tsublist and vsublist respectively. Compute MSE
                // and store it in myContainer
                List<Double> realFreqs_tsublist = new ArrayList<Double>();
                int sum = gts3.size();
                double SE = 0;
                double oneError;
                for (int i=0; i<=sum-1;i++) {
                    realFreqs_tsublist.add(1.0/sum);
                    oneError = realFreqs_tsublist.get(i)-probList.get(i);
                    SE += oneError*oneError;
                }
                myContainer_t.add(SE/sum); // MSE = SE/sum
                myContainer_t_ML.add(curBestScore);   // Likelihood values

                List<Double> realFreqs_vsublist = new ArrayList<Double>();
                sum = gts2.size();
                SE = 0;
                for (int i = 0; i <= sum-1;i++) {
                    realFreqs_vsublist.add(1.0/sum);
                    oneError = realFreqs_vsublist.get(i)-probList.get(i);
                    SE += oneError*oneError;
                }
                myContainer_v.add(SE/sum);
                double ML_v = 0;
                for(int i=0; i<=probList.size()-1; i++){
                    ML_v += Math.log(probList.get(i));
                }
                myContainer_v_ML.add(ML_v);
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

            // Likelihood counterparts
            double total_ML_t = 0;
            double total_ML_v = 0;
            for (double entry: myContainer_t_ML)
                total_ML_t += entry;
            myContainer_t_ML.clear();

            for (double entry: myContainer_v_ML)
                total_ML_v += entry;
            myContainer_v_ML.clear();

            myContainer2_t_ML.add(total_ML_t);  // store total_ML_t
            myContainer2_v_ML.add(total_ML_v);  // store total_ML_v

            for (double entry: myContainer2_t_ML) {
                System.out.println("In myContainer2_t_ML:" + entry);
            }
            for (double entry: myContainer2_v_ML) {
                System.out.println("In myContainer2_v_ML:" + entry);
            }

            // Judgement criteria.
            // Observation: There is randomness in each round of computation such that
            // the third effective digit can be different. Also when trend changing happens, the difference is
            // also pretty small, sometimes around third effective digit.  Therefore, a strict turning might not
            // be fine. Instead we can add a percentage range to determine the peaking phenomena.

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
                        // System.exit(-1);
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
                //          Jump out of the while loop (break) and skip the next if (set isTurningPointFound to true)

                for (int i=0; i <= myContainer2_v.size()-3;i++) {

                    double MSE1 = myContainer2_v.get(i);
                    double MSE2 = myContainer2_v.get(i+1);
                    double MSE3 = myContainer2_v.get(i+2);
                    if (MSE2 < MSE1 && MSE2 < MSE3)  {// turn point found by definition
                        System.out.println("Turning point found. The right reticulation number is " + (i+1));
                        correctK = i+1;
                        isTurningPointFound = true;
                        break;
                    }
                    else if (Math.abs(MSE2-MSE3)/MSE2 < 1e-2) { // turning point found by sticky area criteria
                        System.out.println("Turning point found within 1% MSE. The right reticulation number is " + (i+1));
                        correctK = i+1;
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
        return correctK;
    }

    private Map<String, String> computeAlleleToSpecies(List<Tree> gts, Map<String,List<String>> species2alleles, Set<String> speciesSet){
        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new LinkedHashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                String species = entry.getKey();
                for(String allele: entry.getValue()){
                    allele2species.put(allele,species);
                }
            }
            speciesSet.addAll(species2alleles.keySet());
        }
        else{
            for(Tree gt: gts){
                for(String leaf: gt.getLeaves()){
                    speciesSet.add(leaf);
                }
            }
        }
        return allele2species;
    }

    private DirectedGraphToGraphAdapter<String,PhyloEdge<String>> getStartNetwork(List<Tree> gts, Map<String,String> allele2species, Network<Object> startingNetwork){
        if(startingNetwork == null){
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
            for(NetNode<Object> parent: node.getParents()){
                node.setParentDistance(parent, Double.NaN);
                node.setParentSupport(parent, Double.NaN);
                node.setParentProbability(parent, Double.NaN);
            }
        }

        String newNetwork = network2String(startingNetwork);

        return makeNetwork(newNetwork);

    }

    private Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }

    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                //System.out.println("Start scoring ...");
                //System.out.println("\n"+network2String(network));
                //long start = System.currentTimeMillis();
                Network<Object> speciesNetwork = networkNew2Old(network);

                double score = findUltrametricOptimalBranchLength(speciesNetwork, distinctTrees, species2alleles);


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
                //System.out.println();
                //System.out.println(network2String(speciesNetwork) + ": "+score);
                //System.out.println();
                //System.out.println("End scoring ..." + (System.currentTimeMillis()-start)/1000.0);
                //System.exit(0);
                return score;
            }
        };
    }

    private double findUltrametricOptimalBranchLength(final Network<Object> speciesNetwork, final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles){
        boolean continueRounds = true; // keep trying to improve network


        Map<NetNode, MutableTuple<List<SpeciesPair>, Integer>> node2constraints = new LinkedHashMap<NetNode, MutableTuple<List<SpeciesPair>, Integer>>();
        Map<SpeciesPair, MutableTuple<List<NetNode>, BitSet>> pairHeight2nodes = new LinkedHashMap<SpeciesPair, MutableTuple<List<NetNode>, BitSet>>();
        computeNodeHeightUpperbound(speciesNetwork, node2constraints, pairHeight2nodes);

        Map<NetNode<Object>, Double> node2height = new LinkedHashMap<NetNode<Object>, Double>();
        Map<NetNode, Integer> node2depth = new LinkedHashMap<NetNode, Integer>();
        Map<NetNode, Integer> node2ID = new LinkedHashMap<NetNode, Integer>();
        int id = 0;
        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            node2ID.put(node, id++);
            if(node.isLeaf()){
                node2height.put(node, 0.0);
                node2depth.put(node, 0);
                continue;
            }
            double upperBound = -1;
            if(node2constraints.containsKey(node)){
                upperBound = node2constraints.get(node)._item1.get(0)._time;
            }
            node2height.put(node, upperBound);
            int maxDepth = 0;
            for(NetNode child: node.getChildren()){
                maxDepth = Math.max(maxDepth, node2depth.get(child));
            }
            node2depth.put(node, maxDepth+1);
        }

        NetNode[] nodeArray = new NetNode[node2ID.size()];
        for(Map.Entry<NetNode, Integer> entry: node2ID.entrySet()){
            nodeArray[entry.getValue()] = entry.getKey();
        }
        boolean[][] M = computeM(speciesNetwork, node2ID);

        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            int nodeID = node2ID.get(node);
            double minParent = Double.MAX_VALUE;
            int minParentDepth = -1;
            double maxChild = Double.MIN_VALUE;
            for(NetNode<Object> relateNode: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
                int relateNodeID = node2ID.get(relateNode);
                if(M[relateNodeID][nodeID]){
                    double parentHeight = node2height.get(nodeArray[relateNodeID]);
                    if(parentHeight>0){
                        if(minParent > parentHeight){
                            minParent = parentHeight;
                            minParentDepth = node2depth.get(nodeArray[relateNodeID]);
                        }
                    }
                }
                else if(M[nodeID][relateNodeID]){
                    double childHeight = node2height.get(nodeArray[relateNodeID]);
                    if(childHeight>=0){
                        maxChild = Math.max(maxChild, childHeight);
                    }
                    else{
                        throw new RuntimeException();
                    }
                }
            }
            double currentHeight = node2height.get(node);
            if(currentHeight >= minParent || (currentHeight==-1 && minParent!=Double.MAX_VALUE)){
                int depthDiff = minParentDepth - node2depth.get(node) + 1;
                currentHeight = maxChild + (minParent - maxChild)/depthDiff;
                node2height.put(node, currentHeight);
            }
            else if(currentHeight==-1 && minParent==Double.MAX_VALUE){
                currentHeight = maxChild + 1;
                node2height.put(node, currentHeight);
            }

        }

        double overallMin = Double.MAX_VALUE;
        for(double height: node2height.values()){
            if(height > 0){
                overallMin = Math.min(overallMin, height);
            }
        }
        overallMin = Math.min(overallMin/5, Double.MIN_VALUE);

        //overallMin = 0;

        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            if(node.isLeaf())continue;
            double updatedHeight = node2height.get(node) - overallMin;
            node2height.put(node, updatedHeight);
            for(NetNode child: node.getChildren()){
                child.setParentDistance(node, updatedHeight - node2height.get(child));
                if(child.isNetworkNode()){
                    child.setParentProbability(node,0.5);
                }
            }
        }

        for(NetNode<Object> node: speciesNetwork.bfs()){
            double height = node2height.get(node);
            if(height<0){
                throw new RuntimeException();
            }
            for(NetNode child: node.getChildren()){
                if(height <= node2height.get(child)){
                    throw new RuntimeException();
                }
            }
        }
        //System.out.println(network2String(speciesNetwork));

        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(computeProbability(speciesNetwork, distinctTrees, species2alleles));  // records the GTProb of the network at all times
        final Container<Map<NetNode<Object>, Double>> node2heightContainer = new Container<Map<NetNode<Object>, Double>>(node2height);


        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();

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

                        double lnProb = computeProbability(speciesNetwork, distinctTrees, species2alleles);

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

            }

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

                //TODO
                MutableTuple<List<SpeciesPair>, Integer> spTuple = node2constraints.get(node);
                int currentIndex = -1;
                int maxIndex = -1;
                boolean hasUpperBound = false;
                boolean canLower = false;
                if(spTuple != null){
                    currentIndex = spTuple._item2;
                    canLower = currentIndex!=0;
                    maxIndex = spTuple._item2;
                    SpeciesPair sp;
                    boolean canHigher = true;
                    while(canHigher && spTuple._item1.size()>maxIndex){
                        sp = spTuple._item1.get(maxIndex);
                        canHigher = pairHeight2nodes.get(sp)._item2.cardinality()>1;
                        if(canHigher){
                            maxIndex++;
                        }
                    }
                    if(spTuple._item1.size()>maxIndex){
                        hasUpperBound = true;
                    }
                }
                double minParent = Double.MAX_VALUE;
                if(!node.isRoot()){
                    for(NetNode<Object> parent: node.getParents()){
                        double parentHeight = node2heightContainer.getContents().get(parent);
                        minParent = Math.min(minParent, parentHeight);
                    }
                }
                else{
                    minParent = minHeight.getContents() + _maxBranchLength;
                }

                if(hasUpperBound){
                    maxHeight.setContents(Math.min(minParent, spTuple._item1.get(maxIndex)._time));
                }
                else{
                    maxHeight.setContents(minParent);
                }
                if(canLower){
                    canLower = spTuple._item1.get(currentIndex-1)._time>minHeight.getContents();
                }

                //System.out.println("\nChanging node " + node.getName() + " from " + minHeight.getContents() + " to " + maxHeight.getContents());
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


                        double lnProb = computeProbability(speciesNetwork, distinctTrees, species2alleles);
                        //System.out.print("suggest: "+ suggestedHeight + " " + lnProb + " ");
                        if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                        {
                            lnGtProbOfSpeciesNetwork.setContents(lnProb);
                            node2heightContainer.getContents().put(node, suggestedHeight);
                            //System.out.println( " better ");

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
                            //System.out.println( " worse ");
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

                if(currentIndex!=maxIndex || canLower){
                    double updatedHeight = node2height.get(node);
                    int updatedIndex = 0;
                    for(; updatedIndex<spTuple._item1.size(); updatedIndex++){
                        if(updatedHeight < spTuple._item1.get(updatedIndex)._time){
                            break;
                        }
                    }
                    if(updatedIndex!=currentIndex){
                        if(updatedIndex>currentIndex){
                            for(int i=currentIndex; i<updatedIndex; i++){
                                SpeciesPair sp = spTuple._item1.get(i);
                                MutableTuple<List<NetNode>, BitSet> changedSP = pairHeight2nodes.get(sp);
                                int offBit = changedSP._item1.indexOf(node);
                                changedSP._item2.set(offBit, false);
                            }
                        }
                        else if(updatedIndex<currentIndex){
                            currentIndex = Math.min(currentIndex, spTuple._item1.size()-1);
                            for(int i=currentIndex; i>=updatedIndex; i--){
                                SpeciesPair sp = spTuple._item1.get(i);
                                MutableTuple<List<NetNode>, BitSet> changedSP = pairHeight2nodes.get(sp);
                                int offBit = changedSP._item1.indexOf(node);
                                changedSP._item2.set(offBit, true);
                            }
                        }
                        spTuple._item2 = updatedIndex;
                    }
                }
                //System.out.println(network2String(speciesNetwork) + " " + lnGtProbOfSpeciesNetwork.getContents());
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

        //System.out.println(computeProbability(speciesNetwork, distinctTrees, species2alleles));
        /*
        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(speciesNetwork)){
            if(node.isLeaf()){
                if(species2alleles == null || species2alleles.get(node.getName()).size()<2){
                    node.setParentDistance(node.getParents().iterator().next(), Double.NaN) ;
                }
            }
        }
        */
        //System.out.println("\n"+network2String(speciesNetwork));

        //System.out.println(computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList) + " vs. " + lnGtProbOfSpeciesNetwork.getContents());
        return lnGtProbOfSpeciesNetwork.getContents();
    }

    private void computeNodeHeightUpperbound(Network network, Map<NetNode, MutableTuple<List<SpeciesPair>, Integer>> node2constraints, Map<SpeciesPair, MutableTuple<List<NetNode>, BitSet>> pairHeight2nodes){
        Map<NetNode, Set<String>> node2leaves = new LinkedHashMap<NetNode, Set<String>>();
        for(Object nodeObject: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeObject;
            Set<String> leafSet = new LinkedHashSet<String>();
            if(node.isLeaf()){
                leafSet.add(node.getName());
            }
            else{
                Iterator<NetNode> childIt = node.getChildren().iterator();
                if(node.isNetworkNode()){
                    leafSet.addAll(node2leaves.get(childIt.next()));
                }
                else{
                    NetNode child1 = childIt.next();
                    Set<String> child1Leaves = new LinkedHashSet<String>(node2leaves.get(child1));
                    leafSet.addAll(child1Leaves);
                    NetNode child2 = childIt.next();
                    Set<String> child2Leaves = new LinkedHashSet<String>(node2leaves.get(child2));
                    leafSet.addAll(child2Leaves);
                    Set<String> temp = new LinkedHashSet<String>(child1Leaves);
                    child1Leaves.removeAll(child2Leaves);
                    child2Leaves.removeAll(temp);
                    if(child1Leaves.size()!=0 && child2Leaves.size()!=0){
                        List<SpeciesPair> spList = new ArrayList<SpeciesPair>();
                        for(String leaf1: child1Leaves){
                            for(String leaf2: child2Leaves){
                                SpeciesPair sp = new SpeciesPair(leaf1, leaf2);
                                double time = _pair2time.get(sp);
                                sp.setTime(time);
                                int index = 0;
                                for(SpeciesPair exsp: spList){
                                    if(time < exsp._time){
                                        break;
                                    }
                                    index++;
                                }
                                spList.add(index, sp);
                                MutableTuple<List<NetNode>, BitSet> MutableTuple = pairHeight2nodes.get(sp);
                                if(MutableTuple == null){
                                    MutableTuple = new MutableTuple<List<NetNode>, BitSet>(new ArrayList<NetNode>());
                                    pairHeight2nodes.put(sp, MutableTuple);
                                }
                                MutableTuple._item1.add(node);
                            }
                        }
                        MutableTuple<List<SpeciesPair>, Integer> MutableTuple = new MutableTuple<List<SpeciesPair>, Integer>(spList, 0);
                        node2constraints.put(node, MutableTuple);
                    }
                }
            }
            node2leaves.put(node, leafSet);
        }
        for(MutableTuple<List<NetNode>, BitSet> MutableTuple: pairHeight2nodes.values()){
            BitSet bs = new BitSet(MutableTuple._item1.size());
            bs.set(0, MutableTuple._item1.size(), true);
            MutableTuple.setItem2(bs);
        }
    }

    private boolean[][] computeM(Network<Object> net, Map<NetNode, Integer> node2ID){
        int numNodes = node2ID.size();
        boolean[][] M = new boolean[numNodes][numNodes];
        for(NetNode<Object> node: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(net)){
            int pID = node2ID.get(node);
            //M[pID][pID] = false;
            for(NetNode child: node.getChildren()){
                int cID = node2ID.get(child);
                M[pID][cID] = true;
                for(int i=0; i<numNodes; i++){
                    if(M[cID][i]){
                        M[pID][i] = true;
                    }
                }
            }
        }
        return M;
    }

    private void computePairwiseCoalesceTime(List<Tree> trees, Map<String,String> allele2species, String[] snTaxa){
        for(Tree tr: trees){
            for(TNode n : tr.postTraverse()){
                if(n.isLeaf()){
                    ((STINode<Double>)n).setData(0.0);
                }
                else{
                    double maxTime = 0;
                    for(TNode n_ch : n.getChildren()){
                        maxTime = Math.max(maxTime, ((STINode<Double>) n_ch).getData() + n_ch.getParentDistance());

                        //break;
                    }
                    ((STINode<Double>)n).setData(maxTime);
                }
            }
        }

        Map<STITreeCluster, Double> cluster2height = new LinkedHashMap<STITreeCluster, Double>();

        //Add the cluster containing all taxa to the end of the list.
        STITreeCluster all = new STITreeCluster(snTaxa);
        for (String t : snTaxa) {
            all.addLeaf(t);
        }
        double minTime = -1;
        for(Tree tr:trees){
            TNode root = tr.getRoot();
            if(((STINode<Double>)root).getData() < minTime || minTime == -1){
                minTime = ((STINode<Double>)root).getData();
            }
        }
        cluster2height.put(all, minTime);


        for (Tree tr : trees) {
            for (STITreeClusterWD<Double> tc : tr.getClustersWD(null, true)) {
                STITreeCluster stCluster = new STITreeCluster(snTaxa);
                for (String s : tc.getClusterLeaves()) {
                    String leaf = allele2species==null? s:allele2species.get(s);
                    stCluster.addLeaf(leaf);
                }
                Double preHeight = cluster2height.get(stCluster);
                if(preHeight == null){
                    if(stCluster.getClusterSize()>1){
                        cluster2height.put(stCluster, tc.getData());
                    }
                }
                else{
                    if(preHeight > tc.getData()){
                        cluster2height.put(stCluster, tc.getData());
                    }
                }
            }
        }

        List<STITreeClusterWD<Double>> clusters = new LinkedList<STITreeClusterWD<Double>>();
        for(Map.Entry<STITreeCluster, Double> entry: cluster2height.entrySet()){
            STITreeClusterWD<Double> newCluster = new STITreeClusterWD<Double>(entry.getKey());
            newCluster.setData(entry.getValue());
            int size = entry.getKey().getClusterSize();
            int index = 0;
            for(STITreeClusterWD<Double> cl: clusters){
                if(size > cl.getClusterSize()){
                    break;
                }
                index++;
            }
            clusters.add(index, newCluster);
        }

        //adjust the time
        for(int i=0;i<clusters.size();i++){
            STITreeClusterWD<Double> cl1 = clusters.get(i);
            for(int j=i+1;j<clusters.size();j++){
                STITreeClusterWD<Double> cl2 = clusters.get(j);
                if(cl1.containsCluster(cl2)){
                    if(cl1.getData() < cl2.getData()){
                        cl2.setData(cl1.getData());
                    }
                }
            }
        }

        _pair2time = new LinkedHashMap<SpeciesPair, Double>();
        //ArrayList<STITreeClusterWD<Double>> taxonPairs = new ArrayList<STITreeClusterWD<Double>>();

        for(int i=0;i<snTaxa.length;i++){
            for(int j=i+1;j<snTaxa.length;j++){
                SpeciesPair sp = new SpeciesPair(snTaxa[i], snTaxa[j]);
                _pair2time.put(sp,-1.0);
            }
        }

        //System.out.println(clusters);

        for(STITreeClusterWD<Double> cl1: clusters){
            for(Map.Entry<SpeciesPair, Double> entry: _pair2time.entrySet()){
                if(cl1.containsLeaf(entry.getKey()._species1) && cl1.containsLeaf(entry.getKey()._species2)){
                    if(cl1.getData()<entry.getValue() || entry.getValue()<0){
                        _pair2time.put(entry.getKey(), cl1.getData());
                    }
                }
            }
        }

    }

    private double computeProbability(Network speciesNetwork, List<Tree> geneTrees, Map<String, List<String>> species2alleles) {
        GeneTreeWithBranchLengthProbabilityYF gtp = new GeneTreeWithBranchLengthProbabilityYF();
        List<Double> probList = gtp.calculateGTDistribution(speciesNetwork, geneTrees, species2alleles);
        if(probList.size() < geneTrees.size()){
            throw new RuntimeException("Infinity");
            //return Double.NEGATIVE_INFINITY;
        }
        double total = 0;
        for(Double prob: probList){
            if(prob.isNaN()){
                throw new RuntimeException("Infinity");
            }
            total += Math.log(prob);
        }
        return total;
    }


    private String network2String(Network speciesNetwork){
        RnNewickPrinter<Double> rnNewickPrinter = new RnNewickPrinter<Double>();
        StringWriter sw = new StringWriter();
        rnNewickPrinter.print(speciesNetwork, sw);
        return sw.toString();
    }

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

    class SpeciesPair{
        String _species1;
        String _species2;
        double _time;

        public SpeciesPair(String s1, String s2){
            _species1 = s1;
            _species2 = s2;
        }

        public SpeciesPair(String s1, String s2, double t){
            _species1 = s1;
            _species2 = s2;
            _time = t;
        }

        public void setTime(double t){
            _time = t;
        }

        public int hashCode(){
            return _species1.hashCode() * _species2.hashCode();
        }

        public boolean equals(Object o) {
            if(!(o instanceof SpeciesPair)){
                return false;
            }
            SpeciesPair sp = (SpeciesPair)o;
            if((this._species1.equals(sp._species1) && this._species2.equals(sp._species2)) || (this._species1.equals(sp._species2) && this._species2.equals(sp._species1))){
                return true;
            }
            else{
                return false;
            }
        }

        public String toString(){
            return _species1+"|"+_species2 + ":" + _time;
        }
    }

    class MutableTuple<T1, T2>{
        public T1 _item1;
        public T2 _item2;

        public MutableTuple(T1 item1, T2 item2)
        {
            _item1 = item1;
            _item2 = item2;
        }

        public MutableTuple(T1 item1)
        {
            _item1 = item1;
        }

        public void setItem2(T2 item2){
            _item2 = item2;
        }

        @Override
        public String toString()
        {
            return "(" + _item1 + ", " + _item2 + ")";
        }
    }

}
