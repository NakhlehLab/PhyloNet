package edu.rice.cs.bioinfo.programs.phylonet.algos.treeAugment;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferNetworkPseudoMLFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.HillClimberBase;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.SimpleHillClimbing;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing.SimulatedAnnealingTA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomTopologyNeighborGeneratorTA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.File;
import java.util.*;

/*
 * @author: Zhen Cao
 * @Date: 1/19/2019
 * This class is to infer a species network from gene trees under MPL based on a backbone tree.
 * */

public class treeAugmentInferNetworkMPL extends InferNetworkPseudoMLFromGTT_SingleTreePerLocus {
    protected boolean test = true;
    private  int[] _acceptCount = new int[7];
    private int _mode = 0;//0: mpl, 1: aic, 2: bic
    private int _geneTreeNum;

    public treeAugmentInferNetworkMPL(File intermediate){
        _intermediateResultFile = intermediate;
    }
    public treeAugmentInferNetworkMPL(){
        super();
    }

    /**
     *
     * This function is to infer a species network from input data
     *
     */
    public int[] getAcceptCount()
    {
        return _acceptCount;
    }

    /**
     *
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
        _geneTreeNum = originalData.size();
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
        if(maxReticulations == 0){
            _topologyVsParameterOperation[0] = 0;
            _topologyVsParameterOperation[1] = 1;
        }
        NetworkRandomParameterNeighbourGenerator parameterGenerator = getNetworkRandomParameterNeighbourGenerator(dataForNetworkInference, allele2species, singleAlleleSpecies);
        NetworkRandomNeighbourGenerator allNeighboursStrategy = new NetworkRandomNeighbourGenerator(new NetworkRandomTopologyNeighborGeneratorTA(speciesNetwork, _topologyOperationWeights, maxReticulations, _moveDiameter, _reticulationDiameter, _fixedHybrid, _seed), _topologyVsParameterOperation[0], parameterGenerator, _topologyVsParameterOperation[1], _seed);
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
        if(postOptimization) {
            List summarizedGTs = new ArrayList();
            List gtCorrespondence = new ArrayList();
            _fullLikelihoodCalculator.summarizeData(originalData, allele2species, summarizedGTs, gtCorrespondence);
            LinkedList<Tuple<Network, Double>> updatedResult = optimizeResultingNetworks(_fullLikelihoodCalculator, summarizedGTs, gtCorrespondence, species2alleles, singleAlleleSpecies, resultList);
            resultList.clear();
            resultList.addAll(updatedResult);
            _acceptCount = searcher.GetAcceptCount();

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
     *
     * This function is to infer a species network from input data
     *
     * @param originalData          a collection of input data, gene trees or sequence alignments
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param maxReticulations      the maximal number of reticulations in the inferred species network
     * @param numSol                number of solutions requested to return
     * @param postOptimization      whether optimizing branch lengths and inheritance probabilities of the inferred species networks is needed
     * @param resultList            resulting species networks along with their log likelihood
     */
    public void inferNetwork(List originalData, Map<String,List<String>> species2alleles, int maxReticulations, int numSol, boolean postOptimization, LinkedList<Tuple<Network,Double>> resultList, int mode){
        _mode = mode;
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
        if(maxReticulations == 0){
            _topologyVsParameterOperation[0] = 0;
            _topologyVsParameterOperation[1] = 1;
        }
        NetworkRandomParameterNeighbourGenerator parameterGenerator = getNetworkRandomParameterNeighbourGenerator(dataForNetworkInference, allele2species, singleAlleleSpecies);
        NetworkRandomNeighbourGenerator allNeighboursStrategy = new NetworkRandomNeighbourGenerator(new NetworkRandomTopologyNeighborGeneratorTA(speciesNetwork, _topologyOperationWeights, maxReticulations, _moveDiameter, _reticulationDiameter, _fixedHybrid, _seed), _topologyVsParameterOperation[0], parameterGenerator, _topologyVsParameterOperation[1], _seed);
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

        _acceptCount = searcher.GetAcceptCount();


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
     *
     * This function is to compute the pseudo-likelihood. _mode==1:AIC, _mode==2 BIC
     * @param summarizedData        summarized input data
     * @param dataCorrespondences   relationships between the original data and the data in summarizedData
     * @param species2alleles       mapping from species to alleles sampled from it
     * @param singleAlleleSpecies   species that have only one allele sampled from each
     *
     */
    public Func1<Network, Double> getScoreFunction(final List summarizedData, final List dataCorrespondences, final Map<String, List<String>> species2alleles, final Set<String> singleAlleleSpecies){
        return new Func1<Network, Double>() {
            public Double execute(Network speciesNetwork) {
                Double tmp = computeLikelihood(speciesNetwork, species2alleles, summarizedData, dataCorrespondences, singleAlleleSpecies);
                if(_mode == 1 || _mode == 2){
                    int k = speciesNetwork.getEdgeCount();
                    if(_mode == 1){
                        return 2*tmp-2*k;
                    }
                    else{
                        return 2*tmp-k*Math.log(_geneTreeNum);
                    }

                }
                return tmp;
            }
        };
    }

    protected HillClimberBase getSearchingStrategy(Comparator<Double> comparator, NetworkRandomNeighbourGenerator allNeighboursStrategy){
        return new SimulatedAnnealingTA(comparator, allNeighboursStrategy, _seed);
    }



}
