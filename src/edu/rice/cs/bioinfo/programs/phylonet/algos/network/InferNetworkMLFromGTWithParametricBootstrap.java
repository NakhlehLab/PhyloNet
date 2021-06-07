package edu.rice.cs.bioinfo.programs.phylonet.algos.network;


import edu.rice.cs.bioinfo.library.phylogenetics.parametricbootstrap.network.NetworkParametricBootstrap;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkBootstrap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.File;
import java.util.*;

/**
 * Created with Yun Yu
 * Date: 3/25/13
 * Time: 3:58 PM
 *
 * This class is to use parametric bootstrap to infer species networks under maximum likelihood
 * See "Maximum Likelihood Inference of Reticulate Evolutionary Histories”, Proceedings of the National Academy of Sciences, 2014
 */
public abstract class InferNetworkMLFromGTWithParametricBootstrap{
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
    protected double[] _operationWeights = {0.1,0.1,0.15,0.55,0.15,0.15,2.8};
    protected int _numRuns = 10;
    protected int _numThreads = 1;
    protected Long _seed = null;
    protected Set<String> _fixedHybrid = new HashSet<String>();
    protected boolean _optimizeBL = false;
    Func2<Network, Integer, List> _simulator;
    Func1<List, Tuple<Network,Double>> _inference;



    /**
     * This is the main function for inferring a species network
     *
     * @param originalData          a collection of input data, gene trees or sequence alignments
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param maxReticulations      the maximal number of reticulations in the inferred species network
     * @param numRounds             the number of rounds of bootstrap
     * @param postOptimization      whether optimizing branch lengths and inheritance probabilities of the inferred species networks is needed
     * @param simulatorPath         the path for simulator which generates gene trees down a species network.
     *                              It's the path of ms when branch lengths of gene trees are used; otherwise, none is needed.
     * @param estimator             the method for calculating the support value. Should be either "softwired" or "hardwired"
     */
    public Tuple<Network,Double> inferNetwork(List originalData, Map<String,List<String>> species2alleles, int maxReticulations, int numRounds, boolean postOptimization, String simulatorPath, String estimator){
        _inference = getInferenceMethod(species2alleles, maxReticulations, postOptimization);
        Tuple<Network, Double> inferredNetwork = _inference.execute(originalData);
        _simulator = getSimulator(species2alleles, simulatorPath);
        _inference = getInferenceMethod(species2alleles, maxReticulations, false);
        performBootstrapping(inferredNetwork.Item1, originalData.size(), numRounds, estimator);
        return inferredNetwork;
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
     * @param operationWeights      the weights for different moves for network rearrangement during the search
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
        _operationWeights = operationWeights;
        _numRuns = numRuns;
        _seed = seed;
        _optimizeBL = optimizeBL;
    }


    /**
     * This function is to obtain the method used for inferring species network
     *
     * @param species2alleles    mapping from species to alleles sampled from it
     * @param maxReticulations   the maximal number of networks examined during the search; can be used as one of the termination criterion of the search
     */
    abstract Func1<List, Tuple<Network,Double>> getInferenceMethod(final Map<String,List<String>> species2alleles, final int maxReticulations, final boolean postOptimize);



    /**
     * This function is to obtain the simulator for generating gene trees down a species network
     *
     * @param species2alleles    mapping from species to alleles sampled from it
     * @param MSPath             the path of ms which is only needed when branch lengths of gene trees are used in the inference
     */
    abstract Func2<Network, Integer, List> getSimulator(final Map<String, List<String>> species2alleles, final String MSPath);



    /**
     * This function is to do the bootstrapping
     *
     * @param model         the model network
     * @param sampleSize    the number of gene trees
     * @param numRound      the number of rounds for bootstrapping
     * @param estimator             the method for calculating the support value. Should be either "softwired" or "hardwired"
     */
    private Network performBootstrapping(Network model, Integer sampleSize, int numRound, String estimator){
        List<Network> inferredModels = new ArrayList<Network>();
        for(int i=0; i<numRound; i++){
            List sample = _simulator.execute(model, sampleSize);
            inferredModels.add(_inference.execute(sample).Item1);
        }
        NetworkBootstrap nb = new NetworkBootstrap();
        if(estimator.toLowerCase().equals("softwired")){
            nb.computeSoftwiredNetworkBranchSupport(model, inferredModels);
        }
        else{
            nb.computeHardwiredNetworkBranchSupport(model, inferredModels);
        }
        return model;
    }
}
