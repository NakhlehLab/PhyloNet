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
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.SimpleHillClimbing;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing.SimulatedAnnealingSalterPearL;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomTopologyNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.TwoNetworkRandomPDGGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.math.BigDecimal;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class InferNetworkML {
    /*
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
    protected double[] _topologyOperationWeight;
    protected double[] _topologyVsParameterOperation;
    protected int _numThreads;
    protected int _numRuns;
    protected Long _seed;
    protected File resultFile = null;
    protected File logFile = null;
    protected List<Tuple<Network,Double>> _networksTried = new ArrayList<Tuple<Network, Double>>();
*/
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
    protected double[] _topologyOperationWeight = {0.1,0.1,0.15,0.6,0.15,0.15};
    //protected double[] _topologyOperationWeight = {0.05,0.05,0.1,10,0.1};
    protected double[] _topologyVsParameterOperation = {0.3,0.7};
    //protected double[] _topologyVsParameterOperation = {1,0};
    protected int _numRuns = 5;
    protected int _numThreads = 1;
    protected Long _seed = null;
    protected Set<String> _fixedHybrid = new HashSet<String>();
    protected File _logFile = null;
    protected File _intermediateResultFile = null;
    protected boolean _optimizeBL = false;

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

    public void setStartNetwork(Network startNetwork){
        _startNetwork = startNetwork;
    }

    public void setParallel(int numThreads){
        _numThreads = numThreads;
    }

    public void setNumRuns(int numRuns){
        _numRuns = numRuns;
    }

    public void setSeed(long seed){
        _seed = seed;
    }


    public void setSearchParameter(int maxRounds, int maxTryPerBranch, double improvementThreshold, double maxBranchLength, double Brent1, double Brent2, long maxExaminations, int maxFailure, int moveDiameter, int reticulationDiamter, int parallel, Network startNetwork, Set<String> fixedHybrid, double[] operationWeight, int numRuns, Long seed){
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
        _topologyOperationWeight = new double[operationWeight.length-1];
        _topologyVsParameterOperation = new double[2];
        for(int i=0; i<operationWeight.length; i++){
            if(i!=operationWeight.length-1){
                _topologyOperationWeight[i] = operationWeight[i];
                _topologyVsParameterOperation[0] += _topologyOperationWeight[i];
            }
            else{
                _topologyVsParameterOperation[1] = operationWeight[i];
            }
        }
        _numRuns = numRuns;
        _seed = seed;
    }



    public void inferNetwork(List originalData, Map<String,List<String>> species2alleles, int maxReticulations, int numSol, LinkedList<Tuple<Network,Double>> resultList){
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

        Set<String> singleAlleleSpecies = new HashSet<>();
        if(_optimizeBL){
            _topologyVsParameterOperation[0] = 1;
            _topologyVsParameterOperation[1] = 0;
        }
        else{
            findSingleAlleleSpeciesSet(dataForStartingNetwork, allele2species, singleAlleleSpecies);
        }

        //String startingNetwork = getStartNetwork(dataForStartingNetwork, species2alleles, _fixedHybrid, _startNetwork);
        String startingNetwork = _startNetwork.toString();
        dataForStartingNetwork.clear();

        NetworkRandomNeighbourGenerator allNeighboursStrategy = new NetworkRandomNeighbourGenerator(new NetworkRandomTopologyNeighbourGenerator(_topologyOperationWeight, maxReticulations, _moveDiameter, _reticulationDiameter, _fixedHybrid, _seed), _topologyVsParameterOperation[0], new NetworkRandomParameterNeighbourGenerator(singleAlleleSpecies), _topologyVsParameterOperation[1], _seed);
        Comparator<Double> comparator = getDoubleScoreComparator();
        //SimpleHillClimbing searcher = new SimpleHillClimbing(comparator, allNeighboursStrategy);

        SimulatedAnnealingSalterPearL searcher = new SimulatedAnnealingSalterPearL(comparator, allNeighboursStrategy, _seed);
        searcher.setLogFile(_logFile);
        searcher.setResultFile(_intermediateResultFile);

        Func1<Network, Double> scorer = getScoreFunction(dataForNetworkInference, species2alleles, dataCorrespondence);
        Network speciesNetwork = Networks.readNetwork(startingNetwork);
        searcher.search(speciesNetwork, scorer, numSol, _numRuns, _maxExaminations, _maxFailure, _optimizeBL, resultList); // search starts here

        //
/*
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


    abstract protected String getStartNetwork(List dataForStartingNetwork, Map<String,List<String>> species2alleles, Set<String> hybridSpecies, Network<Object> startingNetwork);


    abstract protected void summarizeData(List originalData, Map<String,String> allele2species, List dataForStartingNetwork, List dataForInferNetwork, List dataCorrespondences);





    protected Comparator<Double> getDoubleScoreComparator(){
        return new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };
    }

    protected Func1<Network, Double> getScoreFunction(final List summarizedData, final Map<String, List<String>> species2alleles, final List dataCorrespondences){
        return new Func1<Network, Double>() {
            public Double execute(Network speciesNetwork) {
                return computeLikelihood(speciesNetwork, species2alleles, summarizedData, dataCorrespondences);
            }
        };
    }


    abstract protected double computeLikelihood(final Network<Object> speciesNetwork, final Map<String, List<String>> species2alleles, final List summarizedData, final List dataCorrespondence);

    abstract protected void findSingleAlleleSpeciesSet(List dataForNetworkInference, Map<String, String> allele2species, Set<String> singleAlleleSpecies);





}