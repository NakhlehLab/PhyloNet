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
 */
public class SimpleHillClimbing extends HillClimberBase{
    protected NetworkNeighbourhoodGenerator _networkGenerator;

    //public SimpleHillClimbing() {}

    public SimpleHillClimbing(Comparator<Double> scoreComparator, NetworkNeighbourhoodGenerator generator) {
        super(scoreComparator);
        _networkGenerator = generator;
    }


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


    protected boolean makeAcceptanceDecision(double currentScore, double newScore){
        return compareTwoScores(newScore, currentScore)>0;
    }

    protected void handleAcceptCase(Network currentNetwork, Ref<Double> currentScore, double newScore){
        currentScore.set(newScore);
        clearConsecutiveFailures();
    }

    protected Network handleRejectCase(Network currentNetwork, Ref<Double> currentScore, double newScore){
        _networkGenerator.undo();
        increaseConsecutiveFailures();
        return currentNetwork;
    }


    protected boolean concludeSearch(){
        return reachMaxConsecutiveFailures() || reachMaxExaminations();
    }



}
