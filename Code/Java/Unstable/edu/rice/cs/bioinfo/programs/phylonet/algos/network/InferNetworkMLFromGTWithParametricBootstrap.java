package edu.rice.cs.bioinfo.programs.phylonet.algos.network;


import edu.rice.cs.bioinfo.library.phylogenetics.parametricbootstrap.network.NetworkParametricBootstrap;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkBootstrap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 3/25/13
 * Time: 3:58 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class InferNetworkMLFromGTWithParametricBootstrap{
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

    Func2<Network, Integer, List> _simulator;
    Func1<List, Tuple<Network,Double>> _inference;



    public Tuple<Network,Double> inferNetwork(List originalData, Map<String,List<String>> species2alleles, int maxReticulations, int numRounds, String simulatorPath, String estimator){
        _inference = getInferenceMethod(species2alleles, maxReticulations);
        Tuple<Network, Double> inferredNetwork = _inference.execute(originalData);
        //System.out.println("\n"+network2String(result.Item1) + " : " + result.Item2);
        _simulator = getSimulator(species2alleles, simulatorPath);

        performBootstrapping(inferredNetwork.Item1, originalData.size(), numRounds, estimator);
        return inferredNetwork;
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

    abstract Func1<List, Tuple<Network,Double>> getInferenceMethod(final Map<String,List<String>> species2alleles, final int maxReticulations);

    abstract Func2<Network, Integer, List> getSimulator(final Map<String, List<String>> species2alleles, final String MSPath);


    private Network performBootstrapping(Network model, Integer sampleSize, int numRound, String estimator){
        List<Network> inferredModels = new ArrayList<>();
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