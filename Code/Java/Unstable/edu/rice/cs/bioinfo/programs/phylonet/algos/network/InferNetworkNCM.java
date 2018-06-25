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
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_Rooted;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.SimpleHillClimbing;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomTopologyNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NonUltrametricNetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.File;
import java.util.*;

/**
 * Created by Hunter Tidwell
 * Date: 1/22/18
 *
 * Based on InferNetworkMP by Yun Yu
 * This class is to infer species networks from gene trees under No Common Mechanism model of Tuffley & Steel
 */
public class InferNetworkNCM{
    protected int _maxFailure = 100;
    protected long _maxExaminations = -1;
    protected int _moveDiameter = -1;
    protected int _reticulationDiameter = -1;
    private Set<String> _fixedHybrid = new HashSet<String>();
    protected Network<Object> _startNetwork = null;

    private int _numProcessors = 1;

    protected double[] _operationWeights = {0.1,0.1,0.15,0.55,0.15,0.15};
    private int _numRuns = 5;
    protected File _logFile = null;
    protected File _intermediateResultFile = null;
    private Long _seed = null;


    /**
     * This function is to set the log file, which saves every network tried along with their likelihood during the search
     */
    public void setLogFile(File file){
        if(!file.exists()){
            try {
                file.createNewFile();
            }catch (Exception e){
                System.err.println(e.getMessage());
                e.getStackTrace();
            }
        }
        _logFile = file;
    }


    /**
     * This function is to set the intermediate result file, which saves the results from every run
     */
    public void setIntermediateResultFile(File file){
        if(!file.exists()){
            try {
                file.createNewFile();
            }catch (Exception e){
                System.err.println(e.getMessage());
                e.getStackTrace();
            }
        }
        _intermediateResultFile = file;
    }

    /**
     * This function is to set all parameters used during the search
     *
     * @param maxExaminations       the maximal number of networks examined during the search; can be used as one of the termination criterion of the search
     * @param maxFailure            the maximal number of consecutive failures during the search before terminating the search; used only in hill climbing
     * @param moveDiameter          the maximal diameter of a move when rearranging a species network
     * @param reticulationDiameter  the maximal diameter of a reticulation in a species network
     * @param numProcessors         the number of threads for parallel computing
     * @param startNetwork          the starting network for search
     * @param fixedHybrid           the species under reticulations in the species network
     * @param operationWeights       the weights for different moves for network rearrangement during the search
     * @param numRuns               the number of independent runs for the search
     * @param seed                  seed for randomness
     */
    public void setSearchParameter(long maxExaminations, int maxFailure, int moveDiameter, int reticulationDiameter, Network startNetwork, Set<String> fixedHybrid, int numProcessors, double[] operationWeights, int numRuns, Long seed){
        _maxExaminations = maxExaminations;
        _moveDiameter = moveDiameter;
        _reticulationDiameter = reticulationDiameter;
        _startNetwork = startNetwork;
        _maxFailure = maxFailure;
        _fixedHybrid = fixedHybrid;
        _numProcessors = 1; //numProcessors; // todo: enable parallelism
        _operationWeights = operationWeights;
        _numRuns = numRuns;
        _seed = seed;
    }


    /**
     * Checks if the user-specified hybrid species are under reticulation nodes in the species network
     */
    private void checkNetworkWithHybrids(Network<Object> startNetwork, Set<String> existingHybrids){
        if(_fixedHybrid.size() == 0) return;
        Map<NetNode, Set<NetNode>> node2children = new HashMap<>();
        Set<NetNode> reticulationNodes = new HashSet<>();
        for(Object o: Networks.postTraversal(startNetwork)){
            NetNode node = (NetNode)o;
            Set<NetNode> childrenSet = new HashSet<>();
            node2children.put(node, childrenSet);
            if(node.isNetworkNode()){
                reticulationNodes.add(node);
            }
            for(Object child: node.getChildren()){
                childrenSet.add((NetNode)child);
                childrenSet.addAll(node2children.get(child));
            }
        }
        for(NetNode reticulation: reticulationNodes){
            for(NetNode child: node2children.get(reticulation)){
                if(child.isLeaf()){
                    if(!_fixedHybrid.contains(child.getName())) {
                        throw new IllegalArgumentException("The starting network contains hybrid that is not in the specified hybrid set.");
                    }else{
                        existingHybrids.add(child.getName());
                    }
                }
            }
        }
    }


    /**
     * This function is to infer a species network from a collection of gene trees
     *
     * @param gts                   the gene trees used to infer the species networks
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param maxReticulations      the maximal number of reticulations in the inferred species network
     * @param numSol                number of solutions requested to return
     * @param blParam               parameter of the prior distribution of branch lengths
     */
    public List<Tuple<Network,Double>> inferNetwork(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles, int maxReticulations, int numSol, double blParam){

        // need to reverse the allele map for integrated probability method
        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new HashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                for(String allele: entry.getValue()){
                    allele2species.put(allele, entry.getKey());
                }
            }
        }


        LinkedList<Tuple<Network, Double>> resultList = new LinkedList<>();
        String startingNetwork = getStartNetwork(gts, species2alleles,_fixedHybrid, _startNetwork);
        NetworkRandomNeighbourGenerator allNeighboursStrategy = new NetworkRandomNeighbourGenerator(new NetworkRandomTopologyNeighbourGenerator(_operationWeights, maxReticulations, _moveDiameter, _reticulationDiameter, _fixedHybrid, _seed), 1, new NonUltrametricNetworkRandomParameterNeighbourGenerator(), 0, _seed);
        Comparator<Double> comparator = getDoubleScoreComparator();
        SimpleHillClimbing searcher = new SimpleHillClimbing(comparator, allNeighboursStrategy);
        searcher.setLogFile(_logFile);
        searcher.setIntermediateResultFile(_intermediateResultFile);

        Func1<Network, Double> scorer = getScoreFunction(gts, allele2species, blParam);

        Network speciesNetwork = Networks.readNetwork(startingNetwork);

        searcher.search(speciesNetwork, scorer, numSol, _numRuns, _maxExaminations, _maxFailure, true, resultList); // search starts here

        //To set inheritance probability // todo: rework after enabling multithreading
        if(_numProcessors != 1){
            for(Tuple<Network, Double> tuple: resultList){
//                MDCOnNetworkYF mdc = new MDCOnNetworkYF();
//                mdc.countExtraCoal(tuple.Item1, gts, species2alleles, new int[gts.size()]);
            }
        }
        return resultList;
    }



    /**
     * Creates reticulation for user-specified hybrid species in a species network
     *
     * @param network   the species network
     * @param hybrid    the user-specified hybrid species
     */
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
        Tuple<NetNode,NetNode> sourceEdge = edgeList.get((int) (Math.random() * numEdges));
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


    /**
     * This function is to obtain starting network for the search
     * MDC on trees is used by default
     *
     * @param gts                 gene trees for inferring the starting network
     * @param species2alleles     mapping from species to alleles sampled from it
     * @param hybridSpecies       species under reticulation nodes in the species network
     * @param startingNetwork     starting network if specified by the users
     *
     * @return  starting network
     */
    private String getStartNetwork(List<MutableTuple<Tree,Double>> gts, Map<String,List<String>> species2alleles, Set<String> hybridSpecies, Network<Object> startingNetwork){
        Set<String> existingHybrids = new HashSet<>();
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
            MDCInference_Rooted mdc = new MDCInference_Rooted();
            Solution sol;
            if(allele2species==null){
                sol = mdc.inferSpeciesTree(gts, false, 1, false, true, -1).get(0);
            }
            else{
                sol = mdc.inferSpeciesTree(gts, allele2species, false, 1, false, true, -1).get(0);
            }
            Tree startingTree= Trees.generateRandomBinaryResolution(sol._st);
            startingNetwork = Networks.readNetwork(startingTree.toString());
        }else{
            checkNetworkWithHybrids(startingNetwork, existingHybrids);
        }

        for(String hybrid: hybridSpecies){
            if(!existingHybrids.contains(hybrid))
                createHybrid(startingNetwork, hybrid);
        }

        Networks.removeAllParameters(startingNetwork);
        Networks.autoLabelNodes(startingNetwork);

        String newNetwork = startingNetwork.toString();
        return newNetwork;

    }


    /**
     * This function is to compare two likelihood scores
     * As a log likelihood, higher is better
     */
    private Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }



    /**
     * This function is to compute the parsimony score (the number of extra lineages) resulting from a species network and a collection of gene trees in parallel
     *
     * @param network           the species network
     * @param gts               gene trees for inferring the starting network
     * @param species2alleles   mapping from species to alleles sampled from it
     * @param xls               the results; has one to one correspondence to gts
     */
    private void computeXLParallel(Network network, List<MutableTuple<Tree,Double>> gts, Map<String, List<String>> species2alleles, int[] xls){
        Thread[] myThreads = new Thread[_numProcessors];

        MDCOnNetworkYF mdc = new MDCOnNetworkYF();
        mdc.setParallel(true);
        mdc.processNetwork(network);

        for(int i=0; i<_numProcessors; i++){
            myThreads[i] = new MyThread(network, gts, species2alleles, xls, mdc);
            myThreads[i].start();
        }

        for(int i=0; i<_numProcessors; i++) {
            try {
                myThreads[i].join();
            } catch (InterruptedException ignore) {
            }
        }

    }


    /**
     * This function is to get the function for calculating likelihood score of a candidate network during the search
     *
     * @param gts               the input gene trees
     * @param allele2species   mapping from allele to species which sampled from it
     */
    private Func1<Network, Double> getScoreFunction(final List<MutableTuple<Tree, Double>> gts, final Map<String, String> allele2species, double blParam){

        List<Tree> gts_stripped = new ArrayList<>();
        for (MutableTuple<Tree, Double> tup: gts) {
            gts_stripped.add(tup.Item1);
        }

        return new Func1<Network, Double>() {
            public Double execute(Network network) {
                //int[] scores = new int[gts.size()];
                if(_numProcessors == 1){
                    // todo: handle parallel

//                    MDCOnNetworkYF scorer = new MDCOnNetworkYF();
//                    scorer.countExtraCoal(network, gts, species2alleles, scores);
                }

                GeneTreeProbabilityIntegrated gtp = new GeneTreeProbabilityIntegrated();
                gtp.setBranchLengthExponentialPrior(blParam, blParam); // todo get both params
                double prob = 0;
                List<Double> probs_l = gtp.calculateGTDistribution(network, gts_stripped, allele2species, false);
                for (int i = 0; i < probs_l.size(); i++) {
                    prob += Math.log(probs_l.get(i)) * gts.get(i).Item2; // multiply log probability by weight factor
//                    System.out.println("adding log of " + probs_l.get(i) + " times " + gts.get(i).Item2);
                }
//                System.out.println("calced prob: " + prob);
                return prob;

//                double total = 0;
//                Iterator<MutableTuple<Tree,Double>> weightIt = gts.iterator();
//                for(int score: scores){
//                    total += score * weightIt.next().Item2;
//                }
//                return total;

            }
        };
    }



    /**
     * This class is for computing the score of a species network in parallel
     */
    private class MyThread extends Thread{
        Network _speciesNetwork;
        List<MutableTuple<Tree,Double>> _geneTrees;
        Map<String, List<String>> _species2alleles;
        int[] _xls;
        MDCOnNetworkYF _mdc;


        public MyThread(Network speciesNetwork, List<MutableTuple<Tree,Double>> geneTrees, Map<String, List<String>> species2alleles, int[] xls, MDCOnNetworkYF mdc){
            _speciesNetwork = speciesNetwork;
            _geneTrees = geneTrees;
            _species2alleles = species2alleles;
            _xls = xls;
            _mdc = mdc;
        }

        public void run() {
            _mdc.countExtraCoal(_speciesNetwork, _geneTrees, _species2alleles, _xls);
        }
    }


}
