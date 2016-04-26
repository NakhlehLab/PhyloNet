package edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkNeighbourhoodRandomWalkGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkNeighbourhoodGenerator;

import java.util.Comparator;

/**
 * Created by yunyu on 10/27/14.
 *
 * This class is a subclass of HillClimberBase
 * It implements the simple hill climbing
 *
 * See "Maximum Likelihood Inference of Reticulate Evolutionary Histories", Proceedings of the National Academy of Sciences, 2014
 */

public class SimpleHillClimbing extends HillClimberBase{
    protected NetworkNeighbourhoodGenerator _networkGenerator;


    /**
     * Constructor / Initializer
     */
    public SimpleHillClimbing(Comparator<Double> scoreComparator, NetworkNeighbourhoodGenerator generator) {
        super(scoreComparator);
        _networkGenerator = generator;
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
        _networkGenerator.mutateNetwork(currentNetwork);
        if(scoreEachTopologyOnce()){
            Double previousResult = getPreviousResult(currentNetwork);
            if(previousResult!=null){
                if(printDetails()){
                    System.out.println("Skip: " + previousResult);
                }
                return previousResult;
            }
        }
        double newScore =  getScore.execute(currentNetwork);
        if(printDetails()){
            System.out.println(newScore + ": " +currentNetwork.toString());
        }

        writeLogFile(currentNetwork, newScore);
        return newScore;
    }


    /**
     * This function is to decide whether the proposed species network should be accepted
     */
    protected boolean makeAcceptanceDecision(double currentScore, double newScore){
        return compareTwoScores(newScore, currentScore) > 0;
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
        clearConsecutiveFailures();
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
        incrementConsecutiveFailures();
        return proposedNetwork;
    }



    /**
     * This function is to decide whether the search should be terminated
     */
    protected boolean concludeSearch(){
        return reachMaxConsecutiveFailures() || reachMaxExaminations();
    }



}
