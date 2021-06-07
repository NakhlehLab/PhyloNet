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


import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.GeneticAlgorithm.GeneticAlgorithmInstance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.HillClimberBase;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.SimpleHillClimbing;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing.SimulatedAnnealingSalterPearL;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomTopologyNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.TwoNetworkRandomPDGGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.util.*;

/**
 * Created by Yun Yu
 * Date: 2/11/13
 * Time: 11:40 AM
 *
 * This is an abstract class for inferring species networks using maximum likelihood
 */
public abstract class InferNetworkML {
    protected NetworkLikelihood _likelihoodCalculator;
    protected boolean _printDetails = false;
    protected int _maxRounds = 100;
    protected int _maxTryPerBranch = 100;
    protected double _improvementThreshold = 0.001;
    protected double _maxBranchLength = 6;
    protected double _Brent1 = 0.01;
    protected double _Brent2 = 0.005;
    protected long _maxExaminations = 1000000;
    protected int _maxFailure = 100;
    protected int _moveDiameter = -1;
    protected int _reticulationDiameter = -1;
    protected Network _startNetwork;
    protected double[] _topologyOperationWeights = {0.1,0.1,0.15,0.55,0.15,0.15};
    protected double[] _topologyVsParameterOperation = {0.3,0.7};
    protected int _numRuns = 10;
    protected int _numThreads = 1;
    protected Long _seed = null;
    protected Set<String> _fixedHybrid = new HashSet<String>();
    protected File _logFile = null;
    protected File _intermediateResultFile = null;
    protected boolean _optimizeBL = false;


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
     * This function is to set the starting network for search
     */
    public void setStartNetwork(Network startNetwork){
        _startNetwork = startNetwork;
    }


    /**
     * This function is to set the number of threads for parallel computing
     */
    public void setParallel(int numThreads){
        _numThreads = numThreads;
    }



    /**
     * This function is to set the number of runs for search
     */
    public void setNumRuns(int numRuns){
        _numRuns = numRuns;
    }



    /**
     * This function is to set seed for randomness
     */
    public void setSeed(long seed){
        _seed = seed;
    }



    /**
     * This function is to set all parameters used during the search
     *
     * @param maxRounds             the maximal rounds when using Brent's method to optimize branch lengths and inheritance probabilities of a network topology from gene trees
     * @param maxTryPerBranch       the maximal number trials of updating one branch length per round when using Brent's method to optimize branch lengths and inheritance probabilities of a network topology from gene trees
     * @param improvementThreshold  the threshold of likelihood improvement between rounds when using Brent's method to optimize branch lengths and inheritance probabilities of a network topology from gene trees
     * @param maxBranchLength       the upper bound of branch lengths when using Brent's method to optimize branch lengths and inheritance probabilities of a network topology from gene trees
     * @param Brent1                rel, which is original stopping criterion of Brent’s algorithm for optimizing branch lengths and inheritance probabilities of a network topology from gene trees
     * @param Brent2                abs, which is original stopping criterion of Brent’s algorithm for optimizing branch lengths and inheritance probabilities of a network topology from gene trees
     * @param maxExaminations       the maximal number of networks examined during the search; can be used as one of the termination criterion of the search
     * @param maxFailure            the maximal number of consecutive failures during the search before terminating the search; used only in hill climbing
     * @param moveDiameter          the maximal diameter of a move when rearranging a species network
     * @param reticulationDiameter  the maximal diameter of a reticulation in a species network
     * @param parallel              the number of threads for parallel computing
     * @param startNetwork          the starting network for search
     * @param fixedHybrid           the species under reticulations in the species network
     * @param operationWeights       the weights for different moves for network rearrangement during the search
     * @param numRuns               the number of independent runs for the search
     * @param optimizeBL            whether the branch lengths and inheritance probabilities of the species network need to be optimized for every network met during the search
     * @param seed                  seed for randomness
     */
    public void setSearchParameter(int maxRounds, int maxTryPerBranch, double improvementThreshold, double maxBranchLength, double Brent1, double Brent2, long maxExaminations, int maxFailure, int moveDiameter, int reticulationDiameter, int parallel, Network startNetwork, Set<String> fixedHybrid, double[] operationWeights, int numRuns, boolean optimizeBL, Long seed){
        _maxRounds = maxRounds;
        _maxTryPerBranch = maxTryPerBranch;
        _improvementThreshold = improvementThreshold;
        _maxBranchLength = maxBranchLength;
        _Brent1 = Brent1;
        _Brent2 = Brent2;
        _maxExaminations = maxExaminations;
        _moveDiameter = moveDiameter;
        _reticulationDiameter = reticulationDiameter;
        _startNetwork = startNetwork;
        _maxFailure = maxFailure;
        _numThreads = parallel;
        _fixedHybrid = fixedHybrid;
        _topologyOperationWeights = new double[operationWeights.length-1];
        _topologyVsParameterOperation = new double[2];
        double totalWeight = 0, totalTopologyWeight = 0;
        for(int i=0; i<operationWeights.length; i++){
            if(i!=operationWeights.length-1){
                _topologyOperationWeights[i] = operationWeights[i];
                _topologyVsParameterOperation[0] += _topologyOperationWeights[i];
                totalTopologyWeight += operationWeights[i];
            }
            else{
                _topologyVsParameterOperation[1] = operationWeights[i];
            }
            totalWeight += operationWeights[i];
        }
        for(int i=0; i<_topologyVsParameterOperation.length; i++){
            _topologyVsParameterOperation[i] /= totalWeight;
        }
        for(int i=0; i<_topologyOperationWeights.length; i++){
            _topologyOperationWeights[i] /= totalTopologyWeight;
        }
        _numRuns = numRuns;
        _optimizeBL = optimizeBL;
        _seed = seed;
    }


    /**
     * This function is to infer a species network from input data
     *
     * @param originalData          a collection of input data, gene trees or sequence alignments
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param maxReticulations      the maximal number of reticulations in the inferred species network
     * @param numSol                number of solutions requested to return
     * @param postOptimization      whether optimizing branch lengths and inheritance probabilities of the inferred species networks is needed
     * @param resultList            resulting species networks along with their log likelihood
     */
    public void inferNetwork(List originalData, Map<String,List<String>> species2alleles, int maxReticulations, int numSol, boolean postOptimization, LinkedList<Tuple<Network,Double>> resultList){
        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new HashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                for(String allele: entry.getValue()){
                    allele2species.put(allele, entry.getKey());
                }
            }
        }

        List dataForNetworkInference = new ArrayList();
        List dataForStartingNetwork = new ArrayList();
        List dataCorrespondence = new ArrayList();
        summarizeData(originalData, allele2species, dataForStartingNetwork, dataForNetworkInference, dataCorrespondence);

        String startingNetwork = getStartNetwork(dataForStartingNetwork, species2alleles, _fixedHybrid, _startNetwork);
        dataForStartingNetwork.clear();
        Network speciesNetwork = Networks.readNetwork(startingNetwork);

        Set<String> singleAlleleSpecies = new HashSet<>();
        findSingleAlleleSpeciesSet(speciesNetwork, species2alleles, singleAlleleSpecies);

        if(_optimizeBL){
            _topologyVsParameterOperation[0] = 1;
            _topologyVsParameterOperation[1] = 0;
        }

        NetworkRandomParameterNeighbourGenerator parameterGenerator = getNetworkRandomParameterNeighbourGenerator(dataForNetworkInference, allele2species, singleAlleleSpecies);
        NetworkRandomNeighbourGenerator allNeighboursStrategy = new NetworkRandomNeighbourGenerator(new NetworkRandomTopologyNeighbourGenerator(_topologyOperationWeights, maxReticulations, _moveDiameter, _reticulationDiameter, _fixedHybrid, _seed), _topologyVsParameterOperation[0], parameterGenerator, _topologyVsParameterOperation[1], _seed);
        Comparator<Double> comparator = getDoubleScoreComparator();

        HillClimberBase searcher;
        if(_optimizeBL){
            searcher = new SimpleHillClimbing(comparator, allNeighboursStrategy);
        } else{
            searcher = getSearchingStrategy(comparator, allNeighboursStrategy);
        }
        searcher.setLogFile(_logFile);
        searcher.setIntermediateResultFile(_intermediateResultFile);
        Func1<Network, Double> scorer = getScoreFunction(dataForNetworkInference, dataCorrespondence, species2alleles, singleAlleleSpecies);

        searcher.search(speciesNetwork, scorer, numSol, _numRuns, _maxExaminations, _maxFailure, _optimizeBL, resultList); // search starts here
        if(postOptimization){
            LinkedList<Tuple<Network,Double>> updatedResult = optimizeResultingNetworks(_likelihoodCalculator, dataForNetworkInference, dataCorrespondence, species2alleles, singleAlleleSpecies, resultList);
            resultList.clear();
            resultList.addAll(updatedResult);
        }
/*
        //this is for searching using genetic algorithm
        NetworkRandomTopologyNeighbourGenerator topologyMutator = new NetworkRandomTopologyNeighbourGenerator(_topologyOperationWeight, maxReticulations, _moveDiameter, _reticulationDiameter, _seed);
        NetworkRandomParameterNeighbourGenerator parameterMutator = new NetworkRandomParameterNeighbourGenerator(singleAlleleSpecies);
        TwoNetworkRandomPDGGenerator recombinationMaker = new TwoNetworkRandomPDGGenerator(_seed);
        Comparator<Double> comparator = getDoubleScoreComparator();
        GeneticAlgorithmInstance searcher = new GeneticAlgorithmInstance(25,comparator,parameterMutator, 0.1, topologyMutator, 0.3, recombinationMaker, 0.1, 5);
        searcher.setLogFile(_logFile);
        Func1<Network, Double> scorer = getScoreFunction(dataForNetworkInference, species2alleles, dataCorrespondence);
        Network speciesNetwork = Networks.readNetwork(startingNetwork);
        searcher.search(speciesNetwork, scorer, numSol, _numRuns, 4000, _optimizeBL, resultList); // search starts here
*/
    }


    /**
     * This function is to optimize branch lengths and inheritance probabilities of inferred species networks
     *
     * @param likelihoodCalculator  the likelihood calculator which can be either full or pseudo likelihood
     * @param summarizedGTs         distinct gene trees
     * @param treeCorrespondence    relationships between distinct gene trees in summarizedGTs and the original gene trees set
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param singleAlleleSpecies   species that have only one allele sampled from it
     * @param resultList            the list of networks to be optimized
     *
     * @return optimized species networks along with their log likelihoods
     */
    protected LinkedList<Tuple<Network,Double>> optimizeResultingNetworks(NetworkLikelihood likelihoodCalculator, List summarizedGTs, List treeCorrespondence, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies, LinkedList<Tuple<Network,Double>> resultList){
        LinkedList<Tuple<Network, Double>> optimizedResults = new LinkedList<>();
        for(Tuple<Network, Double> resultTuple: resultList){
            likelihoodCalculator.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _numThreads);
            double prob = likelihoodCalculator.computeLikelihood(resultTuple.Item1, species2alleles, summarizedGTs, treeCorrespondence, singleAlleleSpecies, true);
            Tuple<Network,Double> newResult = new Tuple<>(resultTuple.Item1, prob);
            int index = 0;
            for(Tuple<Network, Double> updatedTuple: optimizedResults){
                if(updatedTuple.Item2 < prob){
                    break;
                }
                else{
                    index++;
                }
            }
            optimizedResults.add(index, newResult);
        }
        return optimizedResults;
    }



    /**
     * Checks if the user-specified hybrid species are under reticulation nodes in the species network
     */
    protected void checkNetworkWithHybrids(Network<Object> startNetwork, Set<String> existingHybrids){
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


    /**
     * This function is to obtain starting network for the search
     *
     * @param dataForStartingNetwork  data for inferring the starting network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param hybridSpecies           species under reticulation nodes in the species network
     * @param startingNetwork         starting network if specified by the users
     *
     * @return  starting network
     */
    abstract protected String getStartNetwork(List dataForStartingNetwork, Map<String,List<String>> species2alleles, Set<String> hybridSpecies, Network<Object> startingNetwork);


    /**
     * This function is to summarize the input data
     *
     * @param originalData              original input data
     * @param allele2species            mapping from allele to species which it is sampled from
     * @param dataForStartingNetwork    data for inferring the starting network
     * @param dataForInferNetwork       (distinct) data used during the search
     * @param dataCorrespondences       relationships between the original data and the data in dataForInferNetwork
     */
    abstract protected void summarizeData(List originalData, Map<String,String> allele2species, List dataForStartingNetwork, List dataForInferNetwork, List dataCorrespondences);



    /**
     * This function is to compare two likelihood scores
     * For log likelihood, the larger the better
     */
    protected Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }



    /**
     * This function is to get the function for calculating likelihood of a candidate network during the search
     *
     * @param summarizedData        summarized input data
     * @param dataCorrespondences   relationships between the original data and the data in summarizedData
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param singleAlleleSpecies   species that have only one allele sampled from it
     */
    protected Func1<Network, Double> getScoreFunction(final List summarizedData, final List dataCorrespondences, final Map<String, List<String>> species2alleles, final Set<String> singleAlleleSpecies){
        return new Func1<Network, Double>() {
            public Double execute(Network speciesNetwork) {
                return computeLikelihood(speciesNetwork, species2alleles, summarizedData, dataCorrespondences, singleAlleleSpecies);
            }
        };
    }


    /**
     * This function is to compute the likelihood score of a given species network
     *
     * @param speciesNetwork        the species network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param summarizedData        summarized input data
     * @param dataCorrespondences   relationships between the original data and the data in summarizedData
     * @param singleAlleleSpecies   species that have only one allele sampled from each
     */
    protected double computeLikelihood(final Network<Object> speciesNetwork, final Map<String, List<String>> species2alleles, final List summarizedData, final List dataCorrespondences, final Set<String> singleAlleleSpecies){
        _likelihoodCalculator.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _numThreads);
        return _likelihoodCalculator.computeLikelihood(speciesNetwork, species2alleles, summarizedData, dataCorrespondences, singleAlleleSpecies, _optimizeBL);
    }


    /**
     * This function is to find the set of branches whose lengths cannot be estimated so that they can be ignored during the inference
     *
     * @param speciesNetwork        the species network
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param singleAlleleSpecies   species that have only one allele sampled from it
     */
    abstract protected void findSingleAlleleSpeciesSet(Network speciesNetwork, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies);



    /**
     * This function is to obtain the random parameter generator for a species network
     *
     * @param dataForNetworkInference   (distinct) data used for network inference
     * @param allele2species            mapping from allele to the species it is sampled from
     * @param singleAlleleSpecies       species that have only one allele sampled from each
     */
    abstract protected NetworkRandomParameterNeighbourGenerator getNetworkRandomParameterNeighbourGenerator(List dataForNetworkInference, Map<String,String> allele2species, Set<String> singleAlleleSpecies);


    /**
     * This function is to obtain the searching strategy
     */
    abstract protected HillClimberBase getSearchingStrategy(Comparator<Double> comparator, NetworkRandomNeighbourGenerator allNeighboursStrategy);

}
