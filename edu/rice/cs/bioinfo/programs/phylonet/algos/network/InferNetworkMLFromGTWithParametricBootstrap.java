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
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 3/25/13
 * Time: 3:58 PM
 * To change this template use File | Settings | File Templates.
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
    protected File _logFile = null;
    protected File _intermediateResultFile = null;
    protected boolean _optimizeBL = false;

    Func2<Network, Integer, List> _simulator;
    Func1<List, Tuple<Network,Double>> _inference;



    public Tuple<Network,Double> inferNetwork(List originalData, Map<String,List<String>> species2alleles, int maxReticulations, int numRounds, boolean postOptimization, String simulatorPath, String estimator){
        _inference = getInferenceMethod(species2alleles, maxReticulations, postOptimization);
        Tuple<Network, Double> inferredNetwork = _inference.execute(originalData);
        //System.out.println("\n"+network2String(result.Item1) + " : " + result.Item2);
        _simulator = getSimulator(species2alleles, simulatorPath);
        _inference = getInferenceMethod(species2alleles, maxReticulations, false);
        performBootstrapping(inferredNetwork.Item1, originalData.size(), numRounds, estimator);
        return inferredNetwork;
    }


    public void setSearchParameter(int maxRounds, int maxTryPerBranch, double improvementThreshold, double maxBranchLength, double Brent1, double Brent2, long maxExaminations, int maxFailure, int moveDiameter, int reticulationDiamter, int parallel, Network startNetwork, Set<String> fixedHybrid, double[] operationWeights, int numRuns, boolean optimizeBL, Long seed){
        _maxRounds = maxRounds;
        _maxTryPerBranch = maxTryPerBranch;
        _improvementThreshold = improvementThreshold;
        _maxBranchLength = maxBranchLength;
        _Brent1 = Brent1;
        _Brent2 = Brent2;
        _maxExaminations = maxExaminations;
        _moveDiameter = moveDiameter;
        _reticulationDiameter = reticulationDiamter;
        _startNetwork = startNetwork;
        _maxFailure = maxFailure;
        _numThreads = parallel;
        _fixedHybrid = fixedHybrid;
        _operationWeights = operationWeights;
        _numRuns = numRuns;
        _seed = seed;
        _optimizeBL = optimizeBL;
    }

    abstract Func1<List, Tuple<Network,Double>> getInferenceMethod(final Map<String,List<String>> species2alleles, final int maxReticulations, final boolean optimization);

    abstract Func2<Network, Integer, List> getSimulator(final Map<String, List<String>> species2alleles, final String MSPath);


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
