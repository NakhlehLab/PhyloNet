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
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkWholeNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours.AllNeighboursHillClimberFirstBetter;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours.AllNeighboursHillClimberSteepestAscent;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BfsSearch;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
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
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferILSNetworkUsingBLProbabilistically extends MDCOnNetworkYFFromRichNewickJung {
    private Network[] _optimalNetworks;
    private double[] _optimalScores;
    private int _maxRounds;
    private int _maxTryPerBranch;
    private double _improvementThreshold;
    private double _maxBranchLength;
    private double _Brent1;
    private double _Brent2;
    private Long _maxExaminations;
    private Long _maxFailure;
    private int _diameterLimit;
    private Network _startNetwork;
    protected Set<String> _fixedHybrid;
    protected double[] _operationWeight;
    protected int _numRuns;
    protected int _numThreads;
    private int _batchSize = 1;

    private Long _seed;

    private Map<SpeciesPair, Double> _pair2time;


    public InferILSNetworkUsingBLProbabilistically(){
        super(new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER));
    }


    public void setSearchParameter(int maxRounds, int maxTryPerBranch, double improvementThreshold, double maxBranchLength, double Brent1, double Brent2, Long maxExaminations, Long maxFailure, int diameterLimit, int parallel, Network startNetwork, Set<String> fixedHybrid, double[] operationWeight, int numRuns, Long seed){
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
        _numThreads = parallel;
        _fixedHybrid = fixedHybrid;
        _operationWeight = operationWeight;
        _numRuns = numRuns;
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

    public List<Tuple<Network,Double>> inferNetwork(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles, int maxReticulations, int numSol){
        //TODO
        //_batchSize = gts.size()/(_numThreads*2);

        _optimalNetworks = new Network[numSol];
        _optimalScores = new double[numSol];
        Arrays.fill(_optimalScores, Double.NEGATIVE_INFINITY);

        Set<String> speciesSet = new HashSet<String>();
        Map<String, String> allele2species = computeAlleleToSpecies(gts, species2alleles, speciesSet);

        String[] snTaxa = new String[speciesSet.size()];
        int index = 0;
        for(String species: speciesSet){
            snTaxa[index++] = species;
        }
        computePairwiseCoalesceTime(gts, allele2species, snTaxa);

        String startingNetwork = getStartNetwork(gts, allele2species,_fixedHybrid,_startNetwork);

        for(int i=0; i<_numRuns; i++) {
            System.out.println("Run #" + i + ":");

            DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = makeNetwork(startingNetwork);
            //NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
            //AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);
            NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>> allNeighboursStrategy = new NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>>(_operationWeight, makeNode, makeEdge, _seed);
            AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>, Double> searcher = new AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);
            Func1<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, Double> scorer = getScoreFunction(gts, species2alleles);
            Comparator<Double> comparator = getDoubleScoreComparator();
            //DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = getStartNetwork(gts, species2alleles,startNetwork);
            HillClimbResult<DirectedGraphToGraphAdapter<String, PhyloEdge<String>>, Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations, _maxFailure, _diameterLimit); // search starts here
            //HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations,  _diameterLimit); // search starts here
            //HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, maxExaminations, maxReticulations, new Long(100)); // search starts here
        }

        List<Tuple<Network, Double>> resultList = new ArrayList<Tuple<Network, Double>>();
        for(int i=0; i<numSol; i++){
            if(_optimalNetworks[i]!=null)
                resultList.add(new Tuple<Network, Double>(_optimalNetworks[i], _optimalScores[i]));
        }
        //System.out.println("\n #Networks " + result.ExaminationsCount);
        return resultList;
    }

    private Map<String, String> computeAlleleToSpecies(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles, Set<String> speciesSet){
        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new HashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                String species = entry.getKey();
                for(String allele: entry.getValue()){
                    allele2species.put(allele,species);
                }
            }
            speciesSet.addAll(species2alleles.keySet());
        }
        else{
            for(MutableTuple<Tree,Double> gt: gts){
                for(String leaf: gt.Item1.getLeaves()){
                    speciesSet.add(leaf);
                }
            }
        }
        return allele2species;
    }


    private String getStartNetwork(List<MutableTuple<Tree,Double>> gts, Map<String,String> allele2species, Set<String> hybridSpecies, Network<Object> startingNetwork){
        checkNetworkWithHybrids(startingNetwork);

        if(startingNetwork == null){
            MDCInference_DP mdc = new MDCInference_DP();
            Solution sol;
            if(allele2species==null){
                sol = mdc.inferSpeciesTree(gts, false, 1, false, true, -1).get(0);
            }
            else{
                sol = mdc.inferSpeciesTree(gts, allele2species, false, 1, false, true, -1).get(0);
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

        return network2String(startingNetwork);

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


    private Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }


    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final List<MutableTuple<Tree,Double>> distinctTrees, final Map<String, List<String>> species2alleles){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {

                Network<Object> speciesNetwork = networkNew2Old(network);

                double score = findUltrametricOptimalBranchLength(speciesNetwork, distinctTrees, species2alleles);
                System.out.println(network2String(speciesNetwork) + ": "+score);

                if(score > _optimalScores[_optimalNetworks.length-1]){
                    boolean exist = false;
                    for(int i=0; i<_optimalNetworks.length; i++){
                        if(_optimalNetworks[i]==null)break;
                        NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
                        if(metric.computeDistanceBetweenTwoNetworks(speciesNetwork, _optimalNetworks[i]) == 0){
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

                        System.out.println("Optimal Networks");
                        for(int i=0; i<_optimalNetworks.length; i++){
                            if(_optimalNetworks[i]!=null)
                                System.out.println(_optimalScores[i] + ": " + network2String(_optimalNetworks[i]));
                        }
                    }
                }
                System.out.println();
                //System.out.println(network2String(speciesNetwork) + ": "+score);
                //System.out.println();
                //System.out.println("End scoring ..." + (System.currentTimeMillis()-start)/1000.0);
                //System.exit(0);
                return score;
            }
        };
    }


    private void initializeNetwork(Network<Object> speciesNetwork, Map<NetNode, MutableTuple<List<SpeciesPair>, Integer>> node2constraints, Map<NetNode<Object>, Double> node2height){
        Map<NetNode, Integer> node2depth = new Hashtable<NetNode, Integer>();
        Map<NetNode, Integer> node2ID = new Hashtable<NetNode, Integer>();
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
                upperBound = node2constraints.get(node).Item1.get(0)._time;
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

        /*
        double overallMin = Double.MAX_VALUE;
        for(double height: node2height.values()){
            if(height > 0){
                overallMin = Math.min(overallMin, height);
            }
        }
        */
        //overallMin = Math.min(overallMin/5, 0.001);

        double overallMin = 0.0001;

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

    }


    private double findUltrametricOptimalBranchLength(final Network<Object> speciesNetwork, final List<MutableTuple<Tree,Double>> gts, final Map<String, List<String>> species2alleles){
        boolean continueRounds = true; // keep trying to improve network


        Map<NetNode, MutableTuple<List<SpeciesPair>, Integer>> node2constraints = new Hashtable<NetNode, MutableTuple<List<SpeciesPair>, Integer>>();
        Map<SpeciesPair, MutableTuple<List<NetNode>, BitSet>> pairHeight2nodes = new Hashtable<SpeciesPair, MutableTuple<List<NetNode>, BitSet>>();
        computeNodeHeightUpperbound(speciesNetwork, node2constraints, pairHeight2nodes);

        Map<NetNode<Object>, Double> node2height = new Hashtable<NetNode<Object>, Double>();
        initializeNetwork(speciesNetwork, node2constraints, node2height);
        //System.out.println("\n"+network2String(speciesNetwork));
        double initialProb;
        if(_numThreads == 1){
            initialProb = computeProbability(speciesNetwork, gts, species2alleles);
        }else{
            initialProb = computeProbabilityParallel(speciesNetwork, gts, species2alleles);
        }

        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(initialProb);  // records the GTProb of the network at all times
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

                        double lnProb;
                        if(_numThreads == 1){
                            lnProb = computeProbability(speciesNetwork, gts, species2alleles);
                        }else{
                            lnProb = computeProbabilityParallel(speciesNetwork, gts, species2alleles);
                        }

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
                    currentIndex = spTuple.Item2;
                    canLower = currentIndex!=0;
                    maxIndex = spTuple.Item2;
                    SpeciesPair sp;
                    boolean canHigher = true;
                    while(canHigher && spTuple.Item1.size()>maxIndex){
                        sp = spTuple.Item1.get(maxIndex);
                        canHigher = pairHeight2nodes.get(sp).Item2.cardinality()>1;
                        if(canHigher){
                            maxIndex++;
                        }
                    }
                    if(spTuple.Item1.size()>maxIndex){
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
                    maxHeight.setContents(Math.min(minParent, spTuple.Item1.get(maxIndex)._time));
                }
                else{
                    maxHeight.setContents(minParent);
                }
                if(canLower){
                    canLower = spTuple.Item1.get(currentIndex-1)._time>minHeight.getContents();
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


                        double lnProb;
                        if(_numThreads == 1){
                            lnProb = computeProbability(speciesNetwork, gts, species2alleles);
                        }else{
                            lnProb = computeProbabilityParallel(speciesNetwork, gts, species2alleles);
                        }
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
                    for(; updatedIndex<spTuple.Item1.size(); updatedIndex++){
                        if(updatedHeight < spTuple.Item1.get(updatedIndex)._time){
                            break;
                        }
                    }
                    if(updatedIndex!=currentIndex){
                        if(updatedIndex>currentIndex){
                            for(int i=currentIndex; i<updatedIndex; i++){
                                SpeciesPair sp = spTuple.Item1.get(i);
                                MutableTuple<List<NetNode>, BitSet> changedSP = pairHeight2nodes.get(sp);
                                int offBit = changedSP.Item1.indexOf(node);
                                changedSP.Item2.set(offBit, false);
                            }
                        }
                        else if(updatedIndex<currentIndex){
                            currentIndex = Math.min(currentIndex, spTuple.Item1.size()-1);
                            for(int i=currentIndex; i>=updatedIndex; i--){
                                SpeciesPair sp = spTuple.Item1.get(i);
                                MutableTuple<List<NetNode>, BitSet> changedSP = pairHeight2nodes.get(sp);
                                int offBit = changedSP.Item1.indexOf(node);
                                changedSP.Item2.set(offBit, true);
                            }
                        }
                        spTuple.Item2 = updatedIndex;
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
        Map<NetNode, Set<String>> node2leaves = new Hashtable<NetNode, Set<String>>();
        for(Object nodeObject: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeObject;
            Set<String> leafSet = new HashSet<String>();
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
                    Set<String> child1Leaves = new HashSet<String>(node2leaves.get(child1));
                    leafSet.addAll(child1Leaves);
                    NetNode child2 = childIt.next();
                    Set<String> child2Leaves = new HashSet<String>(node2leaves.get(child2));
                    leafSet.addAll(child2Leaves);
                    Set<String> temp = new HashSet<String>(child1Leaves);
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
                                    MutableTuple = new MutableTuple<List<NetNode>, BitSet>(new ArrayList<NetNode>(),null);
                                    pairHeight2nodes.put(sp, MutableTuple);
                                }
                                MutableTuple.Item1.add(node);
                            }
                        }
                        MutableTuple<List<SpeciesPair>, Integer> MutableTuple = new MutableTuple<List<SpeciesPair>, Integer>(spList, 0);
                        node2constraints.put(node, MutableTuple);
                    }
                }
            }
            node2leaves.put(node, leafSet);
        }
        for(MutableTuple<List<NetNode>, BitSet> tuple: pairHeight2nodes.values()){
            BitSet bs = new BitSet(tuple.Item1.size());
            bs.set(0, tuple.Item1.size(), true);
            tuple.Item2 = bs;
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


    private void computePairwiseCoalesceTime(List<MutableTuple<Tree,Double>> trees, Map<String,String> allele2species, String[] snTaxa){
        for(MutableTuple<Tree,Double> tuple: trees){
            Tree tr = tuple.Item1;
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

        Map<STITreeCluster, Double> cluster2height = new Hashtable<STITreeCluster, Double>();

        //Add the cluster containing all taxa to the end of the list.
        STITreeCluster all = new STITreeCluster(snTaxa);
        for (String t : snTaxa) {
            all.addLeaf(t);
        }
        double minTime = -1;
        for(MutableTuple<Tree,Double> tuple: trees){
            Tree tr = tuple.Item1;
            TNode root = tr.getRoot();
            if(((STINode<Double>)root).getData() < minTime || minTime == -1){
                minTime = ((STINode<Double>)root).getData();
            }
        }
        cluster2height.put(all, minTime);


        for(MutableTuple<Tree,Double> tuple: trees){
            Tree tr = tuple.Item1;
            for (STITreeCluster<Double> tc : tr.getClusters(null, true)) {
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

        List<STITreeCluster<Double>> clusters = new LinkedList<STITreeCluster<Double>>();
        for(Map.Entry<STITreeCluster, Double> entry: cluster2height.entrySet()){
            STITreeCluster<Double> newCluster = new STITreeCluster<Double>(entry.getKey());
            newCluster.setData(entry.getValue());
            int size = entry.getKey().getClusterSize();
            int index = 0;
            for(STITreeCluster<Double> cl: clusters){
                if(size > cl.getClusterSize()){
                    break;
                }
                index++;
            }
            clusters.add(index, newCluster);
        }

        //adjust the time
        for(int i=0;i<clusters.size();i++){
            STITreeCluster<Double> cl1 = clusters.get(i);
            for(int j=i+1;j<clusters.size();j++){
                STITreeCluster<Double> cl2 = clusters.get(j);
                if(cl1.containsCluster(cl2)){
                    if(cl1.getData() < cl2.getData()){
                        cl2.setData(cl1.getData());
                    }
                }
            }
        }

        _pair2time = new Hashtable<SpeciesPair, Double>();
        //ArrayList<STITreeCluster<Double>> taxonPairs = new ArrayList<STITreeCluster<Double>>();

        for(int i=0;i<snTaxa.length;i++){
            for(int j=i+1;j<snTaxa.length;j++){
                SpeciesPair sp = new SpeciesPair(snTaxa[i], snTaxa[j]);
                _pair2time.put(sp,-1.0);
            }
        }

        //System.out.println(clusters);

        for(STITreeCluster<Double> cl1: clusters){
            for(Map.Entry<SpeciesPair, Double> entry: _pair2time.entrySet()){
                if(cl1.containsLeaf(entry.getKey()._species1) && cl1.containsLeaf(entry.getKey()._species2)){
                    if(cl1.getData()<entry.getValue() || entry.getValue()<0){
                        _pair2time.put(entry.getKey(), cl1.getData());
                    }
                }
            }
        }

    }



    private class MyThreadFromScratch extends Thread{
        GeneTreeWithBranchLengthProbabilityYF _gtp;
        double[] _probs;


        public MyThreadFromScratch(GeneTreeWithBranchLengthProbabilityYF gtp, double[] probs){
            _probs = probs;
            _gtp = gtp;
        }


        public void run() {
            _gtp.calculateGTDistribution(_batchSize, _probs);

        }
    }

    protected double computeProbability(Network speciesNetwork, List<MutableTuple<Tree,Double>> geneTrees, Map<String, List<String>> species2alleles) {
        GeneTreeWithBranchLengthProbabilityYF gtp = new GeneTreeWithBranchLengthProbabilityYF(speciesNetwork, geneTrees, species2alleles);
        double[] probArray = new double[geneTrees.size()];
        gtp.calculateGTDistribution(_batchSize, probArray);
        double totalProb = 0;
        Iterator<MutableTuple<Tree,Double>> it = geneTrees.iterator();
        for(double prob: probArray){
            totalProb += Math.log(prob) * it.next().Item2;
        }

        return totalProb;
    }

    protected double computeProbabilityParallel(Network<Object> speciesNetwork, List<MutableTuple<Tree,Double>> geneTrees, Map<String, List<String>> species2alleles) {
        double[] probArray = new double[geneTrees.size()];
        Thread[] myThreads = new Thread[_numThreads];

        GeneTreeWithBranchLengthProbabilityYF gtp = new GeneTreeWithBranchLengthProbabilityYF(speciesNetwork, geneTrees, species2alleles);
        gtp.setParallel(true);


        for(int i=0; i<_numThreads; i++){
            myThreads[i] = new MyThreadFromScratch(gtp, probArray);
            myThreads[i].start();
        }

        for(int i=0; i<_numThreads; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        double totalProb = 0;
        Iterator<MutableTuple<Tree,Double>> it = geneTrees.iterator();
        for(double prob: probArray){
            totalProb += Math.log(prob) * it.next().Item2;
        }

        return totalProb;
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

    

}
