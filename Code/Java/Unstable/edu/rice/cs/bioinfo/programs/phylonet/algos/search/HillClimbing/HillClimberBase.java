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
 */
public abstract class HillClimberBase extends SearchBase{
    private boolean _printResultsAfterEachRun = true;

    public HillClimberBase(Comparator<Double> scoreComparator){
        super(scoreComparator);
    }

    private long _examinationsCount = 0;
    protected void incrementExaminations()
    {
        _examinationsCount++;
    }
    protected void clearExaminations()
    {
        _examinationsCount=0;
    }
    protected long getExaminations(){
        return _examinationsCount;
    }
    protected void setExaminations(int count){_examinationsCount = count;}

    private long _maxExaminations = -1;


    protected boolean reachMaxExaminations(){
        if(_maxExaminations != -1){
            return _examinationsCount>=_maxExaminations;
        }
        else{
            return false;
        }
    }

    private int _maxFailures = -1;


    private int _consecutiveFailures;

    protected void clearConsecutiveFailures(){
        _consecutiveFailures = 0;
    }

    protected void increaseConsecutiveFailures(){
        _consecutiveFailures++;
    }

    protected int getConsecutiveFailures(){
        return _consecutiveFailures;
    }

    protected boolean reachMaxConsecutiveFailures(){
        if(_maxFailures != -1){
            return _consecutiveFailures>=_maxFailures;
        }
        else{
            return false;
        }
    }



    public void search(Network startSolution, Func1<Network,Double> getScore, int numOptimums, int numRuns, long maxExaminationsCount, int maxFailures, boolean scoreEachTopologyOnce, LinkedList<Tuple<Network,Double>> resultList){
        _maxExaminations = maxExaminationsCount;
        _maxFailures = maxFailures;
        //_optimalNetworks = resultList;
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
            //System.out.println("Run #" + i);
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

            writeResultFile(i);
        }


    }


    protected void search(Network currentNetwork, Func1<Network,Double> getScore, double initialScore)
    {
        Ref<Double> currentScore = new Ref<>(initialScore);
        while(!concludeSearch())
        {

            if(printDetails()){
                System.out.println("\nTrying #"+getConsecutiveFailures() + " in #" + getExaminations());
            }
            //Ref<Double> newBestScore = new Ref<Double>(null);
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


    protected abstract void handleAcceptCase(Network currentNetwork, Ref<Double> currentScore, double newScore);

    protected abstract Network handleRejectCase(Network currentNetwork, Ref<Double> currentScore, double newScore);

    protected abstract double computeRandomNeighborScore(Network currentNetwork, Func1<Network,Double> getScore);

    protected abstract boolean makeAcceptanceDecision(double originalScore, double newScore);

    protected abstract boolean concludeSearch();
}
