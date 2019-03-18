package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.param.OptimizeAll;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/25/17
 * Time: 1:17 PM
 * To change this template use File | Settings | File Templates.
 *
 * It implements the simulated annealing
 *
 */

//TODO: Move out of MCMC package
public class SimulatedAnnealing {
    public static boolean _printDetails = false;
    private static boolean _enableStrategy = false;
    BiAllelicGTR _BAGTRModel = null;
    private boolean _scoreEachTopologyOnce = false;
    private List<Tuple<Network,Double>> _networksTried = new ArrayList<>();    //only used when _scoreEachTopologyOnce is true
    private LinkedList<Tuple<Network,Double>> _optimalNetworks;
    private LinkedList<Tuple<Network,Double>> _optimalNetworksBack;
    private double _temperature;
    private double _U;
    private double _beta;
    private double _burnin = 10;//1000;
    private double _preScore = -1;
    private boolean _doneBurnin = false;
    private LinkedList<Tuple<Network,MutableTuple<Double,Integer>>> _savedNetworks;
    private int _numSavedNetworks = 10;
    private int _maxUninvestigatedProposals;
    private int _consecutiveUninvestigatedProposals=0;
    private int _investigatedProposalThreshold;
    private boolean _initialized = false;
    private int _numOptimums;
    private int _maxReticulations = 0;
    private int _currentRound = 0;
    private State _state;
    private boolean _printResultsAfterEachRun = true;
    private Random _random = new Random();
    private int _maxRounds = 100;
    protected int _maxTryPerBranch = 100;
    protected double _improvementThreshold = 0.001;
    protected double _maxBranchLength = 6;
    protected double _Brent1 = 0.01;
    protected double _Brent2 = 0.001;
    OptimizeAll _optimizer = null;
    public int count = 0;


    public void setSeed(Long seed) {
        _random = new Random(seed);
    }


    private long _examinationsCount = 0;
    private long _maxExaminations = -1;

    /**
     * This function is to increment the number of networks that have been tried
     */
    protected void incrementExaminations()
    {
        _examinationsCount++;
    }


    /**
     * This function is to reset the number of networks that have been tried
     */
    private void clearExaminations()
    {
        _examinationsCount=0;
    }


    /**
     * This function is to return the number of networks that have been tried
     */
    private long getExaminations(){
        return _examinationsCount;
    }

    /**
     * This function is to check if the number of networks that have been visited has reached the preset maximal
     */
    private boolean reachMaxExaminations(){
        if(_maxExaminations != -1){
            return _examinationsCount >= _maxExaminations;
        }
        else{
            return false;
        }
    }


    private int _maxFailures = -1;
    private int _consecutiveFailures;


    /**
     * This function is to clear the number of consecutive failures
     */
    protected void clearConsecutiveFailures(){
        _consecutiveFailures = 0;
    }


    /**
     * This function is to increment the number of consecutive failures
     */
    protected void incrementConsecutiveFailures(){
        _consecutiveFailures++;
    }


    /**
     * This function is to return the number of consecutive failures
     */
    protected int getConsecutiveFailures(){
        return _consecutiveFailures;
    }



    /**
     * This function is to check if the number of consecutive failures has reached the preset maximal
     */
    protected boolean reachMaxConsecutiveFailures(){
        if(_maxFailures != -1){
            return _consecutiveFailures>=_maxFailures;
        }
        else{
            return false;
        }
    }

    public SimulatedAnnealing() {

    }

    /**
     * This function is to reset the number of consecutive un-investigated proposals
     */
    private void clearConsecutiveUninvestigatedProposals(){
        _consecutiveUninvestigatedProposals = 1;
    }


    /**
     * This function is to increment the number of consecutive un-investigated proposals
     */
    private void incrementConsecutiveUninvestigatedProposals(){
        _consecutiveUninvestigatedProposals++;
    }


    /**
     * This function is to check whether the number of consecutive un-investigated proposals has reached the maximal
     */
    private boolean reachMaxUninvestigatedProposals(){
        return _consecutiveUninvestigatedProposals >= _maxUninvestigatedProposals;
    }

    /**
     * This function is to propose a random neighbor of the current network and compute its score
     *
     * @return  the score of the proposed network
     */
    protected double computeRandomNeighborScore() {
        //if(!_state.getUltrametricNetworkObject().isUltrametric()) {
            //System.out.println("Not ultrametric!!!!!!");
        //}

        _state.propose();
        count++;
        //if(!_state.getUltrametricNetworkObject().isUltrametric()) {
            //System.out.println("Not ultrametric!!!!!!");
        //}

        if(scoreEachTopologyOnce()){
            Double previousResult = getPreviousResult(_state.getNetworkObject());
            if(previousResult!=null){
                if(printDetails()){
                    System.out.println("Skip: " + previousResult);
                }
                return previousResult;
            }
        }
        double newScore = _state.calculateLikelihood();

        boolean needOptimization = _scoreEachTopologyOnce;
        if(needOptimization) {
            _optimizer = new OptimizeAll(_state);
            newScore = _optimizer.propose();
            //_optimizer.propose();
            //_optimizer.undo();
            //_optimizer = null;
        }

        if(printDetails()){
            System.out.println(newScore + ": " +_state.toString());
        }

        /*if(!_doneBurnin){
            if(!Double.isNaN(newScore) && newScore != Utils.INVALID_MOVE){
                _U = Math.max(_U, Math.abs(newScore-_preScore));
                _preScore = newScore;
            }
            if(getExaminations()==_burnin){
                _doneBurnin = true;
                clearExaminations();
                if(printDetails()){
                    System.out.println("Done Burn-in: " + _U);
                    System.out.println();
                }
            }
        }*/
        updateTemperature();
        return newScore;
    }


    /**
     * This function is to update the temperature
     */
    protected void updateTemperature(){
        _temperature = 1.0+10000.0/(1+getExaminations()); //_U/(1+getExaminations()*_beta);////
    }

    /**
     * This is the main function for searching the network space
     *
     */
    public void search(List<MarkerSeq> markerSeqs, BiAllelicGTR BAGTRModel, int numOptimums, int numRuns, long maxExaminationsCount, int maxFailures, boolean scoreEachTopologyOnce, LinkedList<Tuple<Network,Double>> resultList){
        Utils.DISABLE_PARAMETER_MOVES = scoreEachTopologyOnce;
        _maxReticulations = Utils._NET_MAX_RETI;
        _BAGTRModel = BAGTRModel;
        _state = new State(
                Utils._START_NET,
                Utils._START_GT_LIST,
                markerSeqs,
                Utils._POISSON_PARAM,
                Utils._TAXON_MAP,
                BAGTRModel
        );
        _state.calculatePrior();
        _maxExaminations = maxExaminationsCount;
        _maxFailures = maxFailures;
        if(_optimalNetworks!=null){
            resultList.addAll(_optimalNetworks);
        }
        else{
            _optimalNetworks = resultList;
        }
        _optimalNetworksBack = new LinkedList<>();
        _numOptimums = numOptimums;
        _scoreEachTopologyOnce = scoreEachTopologyOnce;
        if(getCashedNetworkListSize()!=0 && _optimalNetworks.size()==0){
            for(Tuple<Network,Double> tried: getAllCashedNetworks()){
                updateOptimalNetworks(_optimalNetworks, tried.Item1, tried.Item2);
            }
        }
        int startingI = 1;
        double score = _state.calculateLikelihood();
        updateOptimalNetworks(_optimalNetworks, _state.getNetworkObject(), score);

        //if(!_state.getUltrametricNetworkObject().isUltrametric()) {
            //System.out.println("Not ultrametric!!!!!!");
        //}

        long tempExaminations = _maxExaminations;
        for(int i=startingI; i<=numRuns; i++){
            _currentRound = i;
            if(i == startingI) {
                Utils._NET_MAX_RETI = _enableStrategy ? 0 : _maxReticulations;
                _maxExaminations = 5000;
            } else {
                _maxExaminations = tempExaminations;
                if(_enableStrategy)
                    Utils._NET_MAX_RETI = Math.min(Utils._NET_MAX_RETI + 1, _maxReticulations);
                else
                    Utils._NET_MAX_RETI = _maxReticulations;
            }
            //_maxExaminations = 10;
            long startingTime = System.currentTimeMillis();
            if(i!=startingI) {
                _examinationsCount = 0;
                _consecutiveFailures = 0;
                incrementExaminations();
            }
            search(score);

            /*if(!scoreEachTopologyOnce){
                LinkedList<Tuple<Network,Double>> updatedResult = optimizeResultingNetworks(resultList);
                resultList.clear();
                resultList.addAll(updatedResult);
            }*/

            if(printDetails() || _printResultsAfterEachRun) {
                System.out.println("\nResults after run #" + i);
                System.out.println("State: " + _state.toString());
                System.out.println("Likelihood: " + _state.calculateLikelihood());
                for (Tuple<Network, Double> result: resultList) {

                    System.out.println(result.Item2 + ": " + Networks.getTopologyString(result.Item1) + " : " + Networks.getFullString(result.Item1));
                }
                System.out.println("Running Time (min): " + (System.currentTimeMillis()-startingTime)/60000.0);
                System.out.println("===============================\n");
            }

            if(false) {
                _state = new State(
                        Utils._START_NET,
                        Utils._START_GT_LIST,
                        markerSeqs,
                        Utils._POISSON_PARAM,
                        Utils._TAXON_MAP,
                        BAGTRModel
                );
            } else {
                Network candidate = _optimalNetworksBack.size() > 0 ? _optimalNetworksBack.get(0).Item1 : null;

                if(candidate == null) {
                    candidate = _optimalNetworks.get(0).Item1;
                    score = _optimalNetworks.get(0).Item2;
                } else {
                    score = _optimalNetworksBack.get(0).Item2;
                }

                _state = new State(
                        Networks.getFullString(candidate),
                        Utils._START_GT_LIST,
                        markerSeqs,
                        Utils._POISSON_PARAM,
                        Utils._TAXON_MAP,
                        BAGTRModel
                );
            }
            _state.calculatePrior();

        }
    }

    /**
     * This function is to search the network space for one round
     *
     */
    protected void search(double initialScore) {
        if (!_initialized) {
            initializeParameters(_state.getNetworkObject(), initialScore);
            _initialized = true;
        }

        //_optimizer = new OptimizeAll(_state);
        //double temp = _optimizer.propose();
        //System.out.println(temp);

        Ref<Double> currentScore = new Ref<>(initialScore);
        while(!concludeSearch())
        {
            if(printDetails()){
                System.out.println("\nTrying #"+getConsecutiveFailures() + " in #" + getExaminations());
            }
            double newScore = computeRandomNeighborScore();
            //System.out.println(_optimalNetworks.get(0).Item2 + ": " + Networks.getTopologyString(_optimalNetworks.get(0).Item1) + " : " + Networks.getFullString(_optimalNetworks.get(0).Item1));
            //if(Networks.getTopologyString(_state.getNetworkObject()).indexOf("(Q)#H") != -1 && Networks.getTopologyString(_state.getNetworkObject()).indexOf("(A)#H") != -1) {
            //if(_state.getNetworkObject().getReticulationCount() >= 2) {
                //System.out.println(_temperature);
                //System.out.println(newScore + "     " + Networks.getTopologyString(_state.getNetworkObject()) + "     " + Networks.getFullString(_state.getNetworkObject()));
            //}

            updateOptimalNetworks(_optimalNetworks, _state.getNetworkObject(), newScore);
            if(_state.getNetworkObject().getReticulationCount() <= Math.min(_currentRound, _maxReticulations) - 1) {
                updateOptimalNetworks(_optimalNetworksBack, _state.getNetworkObject(), newScore);
            }

            updateCashedResults(_state.getNetworkObject(), newScore);
            if (makeAcceptanceDecision(currentScore.get(), newScore)) {
                if (printDetails()) {
                    System.out.println("Accepted ("+currentScore.get()+")\n");
                }
                handleAcceptCase(_state.getNetworkObject(), currentScore, newScore);

            } else {
                if (printDetails()) {
                    System.out.println("Rejected ("+currentScore.get()+")\n");
                }
                handleRejectCase();
            }
            incrementExaminations();
        }
    }


    /**
     * This function is to initialize the temperature
     */
    protected void initializeTemperature(){}



    /**
     * This function is to initialize all parameters
     *
     * @param currentNetwork    the current species network
     * @param initialScore      the score of the initial network
     */
    protected void initializeParameters(Network currentNetwork, double initialScore){
        int ntaxa = currentNetwork.getLeafCount();
        _maxUninvestigatedProposals = 50 * ntaxa;
        _consecutiveUninvestigatedProposals=0;
        _investigatedProposalThreshold = (int) (Math.log(0.05) / (ntaxa * Math.log(1.0 - (1.0 / ntaxa)))) + 1;
        _investigatedProposalThreshold *= ntaxa;
        if (scoreEachTopologyOnce()) {
            _investigatedProposalThreshold *= 2;
        } else {
            _investigatedProposalThreshold *= 8;
        }
        _U = Math.abs(initialScore) / 4;
        setBeta(initialScore, ntaxa, 2000);
        _preScore = initialScore;
        if(_numSavedNetworks!=0){
            _savedNetworks = new LinkedList<>();
        }
        _doneBurnin = false;
        //if(printDetails()){
            System.out.println("Initial setting:");
            System.out.println("MaxUninvestigatedProposals: " + _maxUninvestigatedProposals);
            System.out.println("InvestigatedProposalThreshold: " + _investigatedProposalThreshold);
            System.out.println("U: " + _U);
            System.out.println("Beta: " + _beta);
            System.out.println();
        //}
    }


    /**
     * This function is to decide whether the search should be terminated
     */
    protected boolean concludeSearch(){
        if(printDetails()){
            System.out.println("Consecutive unvestigated network: " + _consecutiveUninvestigatedProposals);
        }

        boolean conclude = reachMaxExaminations() || reachMaxUninvestigatedProposals() || reachMaxConsecutiveFailures();
        if(conclude){
            _initialized = false;
        }
        return conclude;
    }


    /**
     * This function is to set parameter beta
     */
    protected void setBeta(double lnl, int n, int m){

        double c = 10000.0;//0.5;
        double alpha = 0.5;
        _beta = c/((1-alpha)*n + alpha*(-lnl)/m);

        //_beta = 0.05;
    }



    /**
     * This function is to update the saved networks
     *
     * @param newProposal   the proposed species network
     * @param score         the corresponding score of the species network newProposal
     */
    private void updateSavedNetworks(Network newProposal, double score){
        double worstScore = 0;
        int worstIndex = -1;
        int index = 0;
        for(Tuple<Network,MutableTuple<Double,Integer>> network: _savedNetworks){
            if(Networks.hasTheSameTopology(newProposal,network.Item1)){
                if(compareTwoScores(score,network.Item2.Item1)>0){
                    network.Item2.Item1 = score;
                    network.Item2.Item2 = 1;
                    if(printDetails()){
                        System.out.println("Exist in saved networks but higher prob");
                    }
                    clearConsecutiveUninvestigatedProposals();
                    return;
                }
                else{
                    network.Item2.Item2++;
                    if(network.Item2.Item2>_investigatedProposalThreshold){
                        if(printDetails()){
                            System.out.println("Exist in saved networks & exceed the threshold ");
                        }
                        incrementConsecutiveUninvestigatedProposals();
                        return;
                    }
                    else{
                        clearConsecutiveUninvestigatedProposals();
                        if(printDetails()){
                            System.out.println("Exist in saved networks: " + network.Item2.Item1);
                        }
                        return;
                    }
                }
            }
            if(index==0 || compareTwoScores(worstScore, network.Item2.Item1)>0){
                worstIndex = index;
                worstScore = network.Item2.Item1;
            }
            index++;
        }

        if(_savedNetworks.size()==_numSavedNetworks){
            _savedNetworks.remove(worstIndex);
        }
        _savedNetworks.add(new Tuple<Network, MutableTuple<Double, Integer>>(newProposal.clone(), new MutableTuple<Double, Integer>(score,1)));
        if(printDetails()){
            System.out.println("Add to saved networks");
        }
        clearConsecutiveUninvestigatedProposals();
    }


    /**
     * This function is to handle the cases where the proposed species network is accepted
     *
     * @param proposedNetwork   the proposed species network
     * @param currentScore      the score of the current species network
     *                          the proposed species network is a neighbor of the current network
     * @param newScore          the score of the proposed species network
     */
    protected void handleAcceptCase(Network proposedNetwork, Ref<Double> currentScore, double newScore){
        proposedNetwork = proposedNetwork.clone();
        _optimizer = null;
        currentScore.set(newScore);
        if(_numSavedNetworks!=0) {
            updateSavedNetworks(proposedNetwork, newScore);
        }
        //double score = _state.calculateLikelihood();

        if(_scoreEachTopologyOnce && _optimizer != null) {
            _optimizer.undo();
            _optimizer = null;
        }
        _state.accept(Utils.INVALID_MOVE);
        clearConsecutiveFailures();
    }


    /**
     * This function is to handle the cases where the proposed species network is rejected
     *
     */
    private void handleRejectCase(){
        if(_scoreEachTopologyOnce && _optimizer != null) {
            _optimizer.undo();
            _optimizer = null;
        }
        _state.undo(Utils.INVALID_MOVE);
        incrementConsecutiveFailures();
        incrementConsecutiveUninvestigatedProposals();
    }

    /**
     * This function is to set the option of printing details
     */
    protected boolean printDetails(){
        return _printDetails;
    }


    /**
     * This function is to set the option of scoring each topology once.
     * It can be used in MP, or ML when branch lengths and inheritance probabilities are optimized for each topology
     */
    protected boolean scoreEachTopologyOnce(){
        return _scoreEachTopologyOnce;
    }

    /**
     * This function is to get the score of the species network if this network has already been visited before
     */
    protected Double getPreviousResult(Network network){
        for(Tuple<Network,Double> networkTried: _networksTried){
            if(Networks.hasTheSameTopology(networkTried.Item1, network)){
                return networkTried.Item2;
            }
        }
        return null;
    }

    /**
     * This function is to decide whether the proposed species network should be accepted
     */
    protected boolean makeAcceptanceDecision(double currentScore, double newScore){
        /*List<String> names = getTaxaNamesUnderReticulation(_state.getNetworkObject());
        if(CoalescenceHistoriesCounting.getNumTaxaUnderReticulation(_state.getNetworkObject()) > 2){
            return false;
        } else if(CoalescenceHistoriesCounting.getNumTaxaUnderReticulation(_state.getNetworkObject()) != names.size()){
            return false;
        }
        for(String name : names) {
            if(!name.startsWith("HYB")) return false;
        }*/

        if(_examinationsCount % 100 == 0) {
            System.out.println("States examined: " + _examinationsCount);
        }

        if(!_state.getUltrametricNetworkObject().isValid()) {
            return false;
        }

        if(printDetails()){
            System.out.println("Temperature: " + _temperature + " move: " + _state.getOperation().getName() + " current score: " + currentScore + " new score: " + newScore);
            System.out.println(Networks.getTopologyString(_state.getNetworkObject()));
            //System.out.println(_state);
        }
        if(compareTwoScores(newScore, currentScore)>0)
        {
            if(printDetails()) {
                System.out.println("Accept");
            }
            return true;
        }
        else{
            if(_temperature==0){
                if(printDetails()){
                    System.out.println("Reject");
                    System.out.println();
                }
                return false;
            }
            double acceptanceRatio = Math.exp((newScore-currentScore)/_temperature);
            double random = _random.nextDouble();
            if(random < 0.01 && _state.getOperation().getName().equals("Delete-Reticulation")) {
                System.out.println("Accept: Delete-Reticulation " + random);
                return true;
            }

            if(random < acceptanceRatio){
                if(printDetails()){
                    System.out.println("Accept: " + random + " < " + acceptanceRatio);
                }
                return true;
            }
            else {
                if(printDetails()){
                    System.out.println("Reject: " + random + " > " + acceptanceRatio);
                }
                return false;
            }
        }
    }

    /**
     * This function is to update the list of optimal networks if needed
     *
     * @param network   the current network that is being visited
     * @param score     the score of the current network
     */
    protected boolean updateOptimalNetworks(LinkedList<Tuple<Network,Double>> _optimalNetworks, Network network, double score){
        if(_optimalNetworks.isEmpty()){
            _optimalNetworks.add(new Tuple<>(network.clone(), score));
            return true;
        }

        if(compareTwoScores(score, _optimalNetworks.getLast().Item2)>0) {
            for (int i = 0; i < _optimalNetworks.size(); i++) {
                Tuple<Network, Double> entry = _optimalNetworks.get(i);
                if (Networks.hasTheSameTopology(network, entry.Item1)) {
                    if (compareTwoScores(score, entry.Item2)>0){
                        if ((i == 0 || (i != 0 && compareTwoScores(_optimalNetworks.get(i-1).Item2, score)>0)) && (i == _optimalNetworks.size() - 1 || (i != _optimalNetworks.size() - 1 && compareTwoScores(_optimalNetworks.get(i + 1).Item2, score)<0))) {
                            _optimalNetworks.set(i, new Tuple<>(network.clone(), score));
                            return true;
                        }
                        _optimalNetworks.remove(i);
                        break;
                    } else {
                        return false;
                    }
                }
            }
            for (int i = 0; i < _optimalNetworks.size(); i++) {
                if(compareTwoScores(score, _optimalNetworks.get(i).Item2)>0){
                    _optimalNetworks.add(i, new Tuple<>(network.clone(), score));
                    break;
                }

            }

            if(_optimalNetworks.size()>_numOptimums){
                _optimalNetworks.removeLast();
            }
            return true;
        }
        else{
            if(_optimalNetworks.size()<_numOptimums){
                for (int i = 0; i < _optimalNetworks.size(); i++) {
                    Tuple<Network, Double> entry = _optimalNetworks.get(i);
                    if (Networks.hasTheSameTopology(network, entry.Item1)) {
                        return false;
                    }
                }
                _optimalNetworks.add(new Tuple<>(network.clone(), score));
                return true;
            }
            else{
                return false;
            }
        }
    }

    LinkedList<Tuple<Network,Double>> optimizeResultingNetworks(LinkedList<Tuple<Network,Double>> resultList) {

        LinkedList<Tuple<Network,Double>> results = new LinkedList<Tuple<Network,Double>>();

        for(Tuple<Network,Double> entry : resultList) {
            UltrametricNetwork network = new UltrametricNetwork(((BniNetwork)entry.Item1).toStringWithRootPop(), null, null, Utils._TAXON_MAP, _BAGTRModel);
            _optimizer = new OptimizeAll(network);
            double newScore = _optimizer.propose();
            if(results.size() == 0) {
                results.add(new Tuple<>(network.getNetwork(), newScore));
            } else {
                for(int i = 0 ; i < results.size() ; i++) {
                    if(newScore < results.get(i).Item2) {
                        results.add(i, new Tuple<>(network.getNetwork(), newScore));
                        break;
                    } else if(i == results.size() - 1) {
                        results.add(new Tuple<>(network.getNetwork(), newScore));
                        break;
                    }
                }
            }
        }


        return results;
    }

    /**
     * This function is to update the list of networks that have been tried along with the likelihood
     *
     * @param network   the current network that is being visited
     * @param score     the score of the current network
     */
    protected void updateCashedResults(Network network, double score) {
        if (_scoreEachTopologyOnce) {
            _networksTried.add(new Tuple<Network, Double>(network.clone(), score));
        }
    }

    /**
     * This function is to compare two scores of two candidate species networks
     *
     * @return 1 if score1 is better; 0 if they are equal; -1 if score2 is better
     */
    protected int compareTwoScores(double score1, double score2){
        return Double.compare(score1, score2);
    }

    /**
     * This function is to return the number of networks that have been tried
     */
    protected int getCashedNetworkListSize(){
        return _networksTried.size();
    }


    /**
     * This function is to return all the networks that have been tried
     */
    protected List<Tuple<Network,Double>> getAllCashedNetworks(){
        return _networksTried;
    }

}
