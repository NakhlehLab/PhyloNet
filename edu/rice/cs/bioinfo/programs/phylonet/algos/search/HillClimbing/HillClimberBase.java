package edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing;

import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SearchBase;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by yunyu on 10/27/14.
 *
 * This class is a subclass of SearchBase
 * It implements the simple hill climbing
 *
 * See "Maximum Likelihood Inference of Reticulate Evolutionary Histories", Proceedings of the National Academy of Sciences, 2014
 */
public abstract class HillClimberBase extends SearchBase{
    private boolean _printResultsAfterEachRun = true;

    /**
     * Constructor / Initializer
     */
    public HillClimberBase(Comparator<Double> scoreComparator){
        super(scoreComparator);
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
    protected void clearExaminations()
    {
        _examinationsCount=0;
    }


    /**
     * This function is to return the number of networks that have been tried
     */
    protected long getExaminations(){
        return _examinationsCount;
    }


    /**
     * This function is to set the number of networks that have been tried
     */
    protected void setExaminations(int count){_examinationsCount = count;}


    /**
     * This function is to check if the number of networks that have been visited has reached the preset maximal
     */
    protected boolean reachMaxExaminations(){
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


    /**
     * This is the main function for searching the network space
     *
     * @param startSolution         the starting network for search
     * @param getScore              the function that computes the score of a candidate species network
     * @param numOptimums           the number of optimal solutions to return
     * @param numRuns               the number of independent runs
     * @param maxExaminationsCount  the maximal number of networks to try during the search; served as stopping criterion
     * @param maxFailures           the maximal number of consecutive failures; served as stopping criterion
     * @param scoreEachTopologyOnce indicate whether each network topology only need to be evaluated once
     * @param resultList            the resulting networks along with their scores
     */
    public void search(Network startSolution, Func1<Network,Double> getScore, int numOptimums, int numRuns, long maxExaminationsCount, int maxFailures, boolean scoreEachTopologyOnce, LinkedList<Tuple<Network,Double>> resultList){
        _maxExaminations = maxExaminationsCount;
        _maxFailures = maxFailures;
        if(_optimalNetworks!=null){
            resultList.addAll(_optimalNetworks);
        }
        else{
            _optimalNetworks = resultList;
        }
        _numOptimums = numOptimums;
        _scoreEachTopologyOnce = scoreEachTopologyOnce;
        if(getCashedNetworkListSize()!=0 && _optimalNetworks.size()==0){
            for(Tuple<Network,Double> tried: getAllCashedNetworks()){
                updateOptimalNetworks(tried.Item1, tried.Item2);
            }
        }
        int startingI = readResultsFromPreviousRuns() + 1;
        double score = getScore.execute(startSolution);
        updateOptimalNetworks(startSolution, score);
        for(int i=startingI; i<=numRuns; i++){
            long startingTime = System.currentTimeMillis();
            if(i!=startingI) {
                _examinationsCount = 0;
                _consecutiveFailures = 0;
                incrementExaminations();
            }
            search(startSolution.clone(), getScore, score);
            if(printDetails() || _printResultsAfterEachRun) {
                System.out.println("\nResults after run #" + i);
                for (Tuple<Network, Double> result: resultList) {
                    System.out.println(result.Item2 + ": " + result.Item1.toString());
                }
                System.out.println("Running Time (min): " + (System.currentTimeMillis()-startingTime)/60000.0);
                System.out.println("===============================\n");
            }
            writeIntermediateResultFile(i);
        }
    }



    /**
     * This function is to search the network space for one round
     *
     * @param currentNetwork        current species network being examined
     * @param getScore              the function that computes the score of a candidate species network
     * @param initialScore          the score of the currentNetwork
     */
    protected void search(Network currentNetwork, Func1<Network,Double> getScore, double initialScore)
    {
        Ref<Double> currentScore = new Ref<>(initialScore);
        while(!concludeSearch())
        {
            if(printDetails()){
                System.out.println("\nTrying #"+getConsecutiveFailures() + " in #" + getExaminations());
            }
            double newScore = computeRandomNeighborScore(currentNetwork, getScore);

            updateOptimalNetworks(currentNetwork, newScore);
            updateCashedResults(currentNetwork, newScore);
            if (makeAcceptanceDecision(currentScore.get(), newScore)) {
                if (printDetails()) {
                    System.out.println("Accepted ("+currentScore.get()+")\n");
                }
                handleAcceptCase(currentNetwork, currentScore, newScore);

            } else {
                if (printDetails()) {
                    System.out.println("Rejected ("+currentScore.get()+")\n");
                }
                currentNetwork = handleRejectCase(currentNetwork, currentScore, newScore);
            }
            incrementExaminations();
        }
    }


    /**
     * This function is to handle the cases where the proposed species network is accepted
     *
     * @param proposedNetwork   the proposed species network
     * @param currentScore      the score of the current species network
     *                          the proposed species network is a neighbor of the current network
     * @param newScore          the score of the proposed species network
     */
    protected abstract void handleAcceptCase(Network proposedNetwork, Ref<Double> currentScore, double newScore);


    /**
     * This function is to handle the cases where the proposed species network is rejected
     *
     * @param proposedNetwork   the proposed species network
     * @param currentScore      the score of the current species network
     *                          the proposed species network is a neighbor of the current network
     * @param newScore          the score of the proposed species network
     */
    protected abstract Network handleRejectCase(Network proposedNetwork, Ref<Double> currentScore, double newScore);



    /**
     * This function is to propose a random neighbor of the current network and compute its score
     *
     * @param currentNetwork    the current species network which will be rearranged in this function
     * @param getScore          the function that computes the score of a candidate species network
     *
     * @return  the score of the proposed network
     */
    protected abstract double computeRandomNeighborScore(Network currentNetwork, Func1<Network,Double> getScore);


    /**
     * This function is to decide whether the proposed species network should be accepted
     */
    protected abstract boolean makeAcceptanceDecision(double originalScore, double newScore);


    /**
     * This function is to decide whether the search should be terminated
     */
    protected abstract boolean concludeSearch();

    /**
     * Function to change printing behavior
     */
    public void setPrintResultsAfterEachRun(boolean doYouWantToPrint) {
        _printResultsAfterEachRun = doYouWantToPrint;
    }
}
