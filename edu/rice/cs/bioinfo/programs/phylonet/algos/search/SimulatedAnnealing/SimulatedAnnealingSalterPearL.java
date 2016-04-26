package edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.util.Comparator;
import java.util.LinkedList;

/**
 * Created by yunyu on 11/4/14.
 *
 * This class is a subclass of SimulatedAnnealingBase.
 *
 * It implements the simulated annealing in Salter and Pearl 2001, which is the same as the one in STEM
 */

public class SimulatedAnnealingSalterPearL extends SimulatedAnnealingBase{
    private double _U;
    private double _beta;
    private double _burnin = 1000;
    private double _preScore = -1;
    private boolean _doneBurnin = false;
    private LinkedList<Tuple<Network,MutableTuple<Double,Integer>>> _savedNetworks;
    private int _numSavedNetworks = 10;
    private int _maxUninvestigatedProposals;
    private int _consecutiveUninvestigatedProposals=0;
    private int _investigatedProposalThreshold;
    private boolean _initialized = false;


    /**
     * Constructor / Initializer
     */
    public SimulatedAnnealingSalterPearL(Comparator<Double> scoreComparator, NetworkNeighbourhoodGenerator generator, Long seed) {
        super(scoreComparator, generator, seed);

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
     * @param currentNetwork    the current species network which will be rearranged in this function
     * @param getScore          the function that computes the score of a candidate species network
     *
     * @return  the score of the proposed network
     */
    protected double computeRandomNeighborScore(Network currentNetwork, Func1<Network, Double> getScore) {
        double newScore = super.computeRandomNeighborScore(currentNetwork, getScore);
        if(!_doneBurnin){
            if(!Double.isNaN(newScore)){
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
        }
        updateTemperature();
        return newScore;
    }


    /**
     * This function is to update the temperature
     */
    protected void updateTemperature(){
        _temperature = _U/(1+getExaminations()*_beta);
    }



    /**
     * This function is to search the network space for one round
     *
     * @param currentNetwork        current species network being examined
     * @param getScore              the function that computes the score of a candidate species network
     * @param initialScore          the score of the currentNetwork
     */
    protected void search(Network currentNetwork, Func1<Network,Double> getScore, double initialScore) {
        if (!_initialized) {
            initializeParameters(currentNetwork, initialScore);
            _initialized = true;
        }
        super.search(currentNetwork, getScore, initialScore);
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
        _maxUninvestigatedProposals = 100 * ntaxa;
        _consecutiveUninvestigatedProposals=0;
        _investigatedProposalThreshold = (int) (Math.log(0.05) / (ntaxa * Math.log(1.0 - (1.0 / ntaxa)))) + 1;
        _investigatedProposalThreshold *= ntaxa;
        if (scoreEachTopologyOnce()) {
            _investigatedProposalThreshold *= 4;
        } else {
            _investigatedProposalThreshold *= 12;
        }
        _U = Math.abs(initialScore) / 4;
        setBeta(initialScore, ntaxa, 2000);
        _preScore = initialScore;
        if(_numSavedNetworks!=0){
            _savedNetworks = new LinkedList<>();
        }
        _doneBurnin = false;
        if(printDetails()){
            System.out.println("Initial setting:");
            System.out.println("MaxUninvestigatedProposals: " + _maxUninvestigatedProposals);
            System.out.println("InvestigatedProposalThreshold: " + _investigatedProposalThreshold);
            System.out.println("U: " + _U);
            System.out.println("Beta: " + _beta);
            System.out.println();
        }
    }


    /**
     * This function is to decide whether the search should be terminated
     */
    protected boolean concludeSearch(){
        if(printDetails()){
            System.out.println("Consecutive unvestigated network: " + _consecutiveUninvestigatedProposals);
        }

        boolean conclude = reachMaxExaminations() || reachMaxUninvestigatedProposals();
        if(conclude){
            _initialized = false;
        }
        return conclude;
    }


    /**
     * This function is to set parameter beta
     */
    protected void setBeta(double lnl, int n, int m){
        /*
        double c = 0.5;
        double alpha = 0.5;
        _beta = c/((1-alpha)*n + alpha*(-lnl)/m);
        */
        _beta = 0.01;
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
        currentScore.set(newScore);
        if(_numSavedNetworks!=0) {
            updateSavedNetworks(proposedNetwork, newScore);
        }
    }


    /**
     * This function is to handle the cases where the proposed species network is rejected
     *
     * @param proposedNetwork   the proposed species network
     * @param currentScore      the score of the current species network
     *                          the proposed species network is a neighbor of the current network
     * @param newScore          the score of the proposed species network
     */
    protected Network handleRejectCase(Network proposedNetwork, Ref<Double> currentScore, double newScore){
        _networkGenerator.undo();
        incrementConsecutiveUninvestigatedProposals();
        return proposedNetwork;
    }



}
