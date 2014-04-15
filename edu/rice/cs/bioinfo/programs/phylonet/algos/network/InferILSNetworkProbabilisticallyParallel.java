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
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferILSNetworkProbabilisticallyParallel extends MDCOnNetworkYFFromRichNewickJung {
    protected Network[] _optimalNetworks;
    protected double[] _optimalScores;
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
    protected double[] _operationWeight;
    protected int _numThreads;
    protected int _numRuns;
    protected Long _seed;


    public InferILSNetworkProbabilisticallyParallel(){
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
        _optimalNetworks = new Network[numSol];
        _optimalScores = new double[numSol];
        Arrays.fill(_optimalScores, Double.NEGATIVE_INFINITY);

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
        summarizeGeneTrees1(gts, allele2species, distinctNonBinaryTrees);

        String startingNetwork = getStartNetwork(distinctNonBinaryTrees, species2alleles, _fixedHybrid, _startNetwork);
        //System.out.println(startingNetwork);
        List<Tree> distinctBinaryTrees = new ArrayList<Tree>();
        List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList = new ArrayList<Tuple<MutableTuple<Tree,Double>, Set<Integer>>>();
        summarizeGeneTrees2(distinctNonBinaryTrees, allele2species, distinctBinaryTrees, nbTreeAndCountAndBinaryIDList);
        //System.out.println("\nDone reading gene trees:" + distinctBinaryTrees.size());


        for(int i=1; i<=_numRuns; i++){
            DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = makeNetwork(startingNetwork);
            //NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkWholeNeighbourhoodGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
            //AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);

            NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> allNeighboursStrategy = new NetworkNeighbourhoodRandomWalkGenerator<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(_operationWeight, makeNode, makeEdge, _seed);
            AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Double> searcher = new AllNeighboursHillClimberFirstBetter<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>, Double>(allNeighboursStrategy);

            Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> scorer = getScoreFunction(distinctBinaryTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
            Comparator<Double> comparator = getDoubleScoreComparator();
            //DirectedGraphToGraphAdapter<String,PhyloEdge<String>> speciesNetwork = getStartNetwork(gts, species2alleles,startNetwork);
            HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations, _maxFailure, _diameterLimit); // search starts here
            //HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Double> result = searcher.search(speciesNetwork, scorer, comparator, _maxExaminations, maxReticulations,  _diameterLimit); // search starts here
/*
            System.out.println("Results after run #" + i);
            for(int j=0; j<_optimalNetworks.length; j++){
                System.out.println(_optimalScores[j] + ": " + network2String(_optimalNetworks[j]));
            }
            System.out.println("===============================\n\n");
*/
        }
        List<Tuple<Network, Double>> resultList = new ArrayList<Tuple<Network, Double>>();
        for(int i=0; i<numSol; i++){
            resultList.add(new Tuple<Network, Double>(_optimalNetworks[i], _optimalScores[i]));
        }
        return resultList;
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
        originalGTs.clear();
    }


    protected void summarizeGeneTrees2(List<MutableTuple<Tree,Double>> originalGTs, Map<String,String> allele2species, List<Tree> distinctGTs, List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList){
        Map<String,Integer> exp2ID = new HashMap<String, Integer>();
        //int temp=0;
        for(MutableTuple<Tree,Double> gt: originalGTs){
            //System.out.println(temp + " :" + distinctGTs.size());
            //temp++;
            Set<Integer> binaryIDs = new HashSet<Integer>();
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
            nbTreeAndCountAndBinaryIDList.add(new Tuple<MutableTuple<Tree,Double>, Set<Integer>>(gt, binaryIDs));

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


    protected Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double> getScoreFunction(final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Double>() {
            public Double execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                Network<Object> speciesNetwork = networkNew2Old(network);
                double score = findNonUltrametricOptimalBranchLength(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
                boolean optimalChanged = false;
                int originalIndex  = -1;
                if(score > _optimalScores[_optimalNetworks.length-1]){

                    for(int i=0; i<_optimalNetworks.length; i++){
                        if(_optimalNetworks[i]==null)break;

                        if(edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeClusterDistance(speciesNetwork, _optimalNetworks[i])[2]<0.000001 &&
                                edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeTripartitionDistance(speciesNetwork, _optimalNetworks[i])[2]<0.000001){
                            originalIndex = i;
                            break;
                        }
                    }
                    if(originalIndex == -1 || score>_optimalScores[originalIndex]){
                        optimalChanged = true;
                        int index = -1;
                        for(int i=0; i<_optimalScores.length; i++){
                            if(score > _optimalScores[i]){
                                index = i;
                                break;
                            }
                        }
                        int startingMovingPoint = originalIndex == -1? _optimalScores.length-1:originalIndex;
                        for(int i=startingMovingPoint; i>index; i--){
                            _optimalNetworks[i] = _optimalNetworks[i-1];
                            _optimalScores[i] = _optimalScores[i-1];
                        }
                        _optimalScores[index] = score;
                        _optimalNetworks[index] = string2Network(network2String(speciesNetwork));
                    }
                }

                System.gc();

                /*
                System.out.println(score +  ": "+ network2String(speciesNetwork));
                if(optimalChanged){
                    System.out.println("Optimal ones:");
                    for(int i=0; i<_optimalNetworks.length; i++){
                        if(_optimalNetworks[i]!=null){
                            System.out.println(_optimalScores[i]+": " + network2String(_optimalNetworks[i]));
                        }
                    }
                }

                System.out.println();
                */
                return score;
            }
        };
    }


    protected double findUltrametricOptimalBranchLength(final Network<Object> speciesNetwork, final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList){
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
            List<Proc> assigmentActions = new ArrayList<Proc>();

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

                assigmentActions.add(new Proc()
                {
                    public void execute()
                    {
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
                });
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
                });
            }

            if(_seed==null) {
                Collections.shuffle(assigmentActions);
            }

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



    protected double findNonUltrametricOptimalBranchLength(final Network<Object> speciesNetwork, final List<Tree> distinctTrees, final Map<String, List<String>> species2alleles, final List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList){
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
        if(_numThreads > 1){
            initalProb = computeProbabilityParallel(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
        }
        else{
            initalProb = computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);
        }
        //System.out.println(initalProb);
        //System.out.print("\n"+(System.currentTimeMillis()-start));
        //computeProbability(speciesNetwork, distinctTrees, species2alleles, nbTreeAndCountAndBinaryIDList);

        //System.out.println();
        //System.out.println(network2String(speciesNetwork));
        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(initalProb);  // records the GTProb of the network at all times
        //final Container<Integer> callCount = new Container<Integer>(0);


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
                                    //System.out.print(" l/"+child.getName()+" ");
                                    //callCount.setContents(callCount.getContents()+1);
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
                                //callCount.setContents(callCount.getContents()+1);
                                //System.out.print(" p/"+child.getName()+" ");
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
            if(_seed==null) {
                Collections.shuffle(assigmentActions);
            }
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


    protected double computeProbabilityParallel(Network<Object> speciesNetwork, List<Tree> distinctTrees, Map<String, List<String>> species2alleles, List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList) {
        double[] probs = new double[distinctTrees.size()];
        Thread[] myThreads = new Thread[_numThreads];

        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        gtp.setParallel(true);
        gtp.preProcess(speciesNetwork, distinctTrees, true);


        for(int i=0; i<_numThreads; i++){
            myThreads[i] = new MyThreadFromScratch(gtp, speciesNetwork, distinctTrees, species2alleles, probs);
            myThreads[i].start();
        }

        for(int i=0; i<_numThreads; i++){
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
        for(Tuple<MutableTuple<Tree,Double>, Set<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double totalProb = 0;
            for(int id: triple.Item2){
                totalProb += probs[id];
            }
            initialProb += Math.log(totalProb) * triple.Item1.Item2;
        }

        /*
        System.out.println("\n" + initialProb);
        System.out.println(System.currentTimeMillis()-start);
        System.exit(0);
        */

        return initialProb;
    }


    protected double computeProbability(Network speciesNetwork, List<Tree> geneTrees, Map<String, List<String>> species2alleles, List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList) {
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        double[] probList = new double[geneTrees.size()];
        gtp.calculateGTDistribution(speciesNetwork, geneTrees, species2alleles, probList);
        double total = 0;
        for(Tuple<MutableTuple<Tree,Double>, Set<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double totalProb = 0;
            for(int id: triple.Item2){
                totalProb += probList[id];
            }
            total += Math.log(totalProb) * triple.Item1.Item2;
        }
        return total;
    }


    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList, NetNode child, NetNode parent) {
        Set<NetNode> childNodes = new HashSet<NetNode>();
        childNodes.add(child);
        Set<NetNode> parentNodes = new HashSet<NetNode>();
        parentNodes.add(parent);
        if(_numThreads>1){
            return computeProbabilityParallel(speciesNetwork, geneTrees, nbTreeAndCountAndBinaryIDList, childNodes, parentNodes);
        }
        else{
            return computeProbability(speciesNetwork, geneTrees, nbTreeAndCountAndBinaryIDList, childNodes, parentNodes);
        }
    }


    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList, NetNode node, boolean changeBranchLength) {
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

        if(_numThreads>1){
            return computeProbabilityParallel(speciesNetwork, geneTrees, nbTreeAndCountAndBinaryIDList, childNodes, parentNodes);
        }
        else{
            return computeProbability(speciesNetwork, geneTrees, nbTreeAndCountAndBinaryIDList, childNodes, parentNodes);
        }
    }





    public double computeProbability(Network speciesNetwork, List<Tree> geneTrees, final List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList, Set<NetNode> childNodes, Set<NetNode> parentNodes) {
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        double[] probList = new double[geneTrees.size()];
        gtp.calculateGTDistribution(speciesNetwork, geneTrees, childNodes, parentNodes, probList);
        double total = 0;
        for(Tuple<MutableTuple<Tree,Double>, Set<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double totalProb = 0;
            for(int id: triple.Item2){
                totalProb += probList[id];
            }
            total += Math.log(totalProb) * triple.Item1.Item2;
        }
        return total;
    }


    public double computeProbabilityParallel(Network<Object> speciesNetwork, List<Tree> distinctTrees, final List<Tuple<MutableTuple<Tree,Double>, Set<Integer>>> nbTreeAndCountAndBinaryIDList, Set<NetNode> childNodes, Set<NetNode> parentNodes) {
        double[] probs = new double[distinctTrees.size()];
        Thread[] myThreads = new Thread[_numThreads];

        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        gtp.setParallel(true);
        gtp.preProcess(speciesNetwork, distinctTrees, false);

        //System.out.println("\ngts:" +distinctTrees);
        for(int i=0; i<_numThreads; i++){
            myThreads[i] = new MyThreadFromNonScratch(gtp, speciesNetwork, distinctTrees, childNodes, parentNodes, probs);
            myThreads[i].start();
        }

        for(int i=0; i<_numThreads; i++){
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {}
        }

        double initialProb = 0;
        for(Tuple<MutableTuple<Tree,Double>, Set<Integer>> triple: nbTreeAndCountAndBinaryIDList){
            double totalProb = 0;
            for(int id: triple.Item2){
                totalProb += probs[id];
            }
            initialProb += Math.log(totalProb) * triple.Item1.Item2;
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
