package edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.Comparator;
import java.util.LinkedList;

/**
 * Created by yunyu on 11/4/14.
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


    public SimulatedAnnealingSalterPearL(Comparator<Double> scoreComparator, NetworkNeighbourhoodGenerator generator, Long seed) {
        super(scoreComparator, generator, seed);
    }

    private void clearConsecutiveUninvestigatedProposals(){
        _consecutiveUninvestigatedProposals = 1;
    }

    private void increaseConsecutiveUninvestigatedProposals(){
        _consecutiveUninvestigatedProposals++;
    }

    private boolean reachMaxUninvestigatedProposals(){
        return _consecutiveUninvestigatedProposals>=_maxUninvestigatedProposals;
    }

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

    protected void updateTemperature(){

        _temperature = _U/(1+getExaminations()*_beta);

    }


    protected void search(Network currentNetwork, Func1<Network,Double> getScore, double initialScore) {
        if (!_initialized) {
            initializeParameters(currentNetwork, initialScore);
            _initialized = true;
        }
        super.search(currentNetwork, getScore, initialScore);
    }


    protected void initializeTemperature(){}

    protected void initializeParameters(Network currentNetwork, double initialScore){
        int ntaxa = currentNetwork.getLeafCount();
        _maxUninvestigatedProposals = 300 * ntaxa;
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

    protected boolean concludeSearch(){
        if(printDetails()){
            System.out.println("Consecutive unvestigated network: " + _consecutiveUninvestigatedProposals);
        }

        boolean conclude = reachMaxExaminations() || reachMaxUninvestigatedProposals();

        if(conclude){
            /*
            System.out.println("Consecutive unvestigated network: " + _consecutiveUninvestigatedProposals);
            System.out.println("Iterations: " + getExaminations());
            for(Tuple<Network,MutableTuple<Double,Integer>> network: _savedNetworks){
                System.out.println(network.toString());
            }
            */
            _initialized = false;
        }

        return conclude;
    }


    protected void setBeta(double lnl, int n, int m){
        /*
        double c = 0.5;
        double alpha = 0.5;
        _beta = c/((1-alpha)*n + alpha*(-lnl)/m);
        */
        _beta = 0.01;
    }

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
                        increaseConsecutiveUninvestigatedProposals();
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


    protected void handleAcceptCase(Network currentNetwork, Ref<Double> currentScore, double newScore){

        currentScore.set(newScore);
        if(_numSavedNetworks!=0) {
            updateSavedNetworks(currentNetwork, newScore);
        }
    }

    protected Network handleRejectCase(Network currentNetwork, Ref<Double> currentScore, double newScore){
        _networkGenerator.undo();
        increaseConsecutiveUninvestigatedProposals();
        return currentNetwork;
    }



}
