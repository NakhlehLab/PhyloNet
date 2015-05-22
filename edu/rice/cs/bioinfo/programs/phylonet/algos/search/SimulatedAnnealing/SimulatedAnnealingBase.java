package edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.SimpleHillClimbing;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkNeighbourhoodGenerator;

import java.util.Comparator;
import java.util.Random;

/**
 * Created by yunyu on 11/3/14.
 */
public abstract class SimulatedAnnealingBase extends SimpleHillClimbing {
    Random _random;
    protected double _temperature = -1;


    public SimulatedAnnealingBase(Comparator<Double> scoreComparator, NetworkNeighbourhoodGenerator generator, Long seed) {
        super(scoreComparator, generator);
        initializeTemperature();
        if(seed!=null) {
            _random = new Random(seed);
        }else{
            _random = new Random();
        }
    }


    protected double computeRandomNeighborScore(Network currentNetwork, Func1<Network, Double> getScore) {
        double newScore = super.computeRandomNeighborScore(currentNetwork, getScore);
        updateTemperature();
        return newScore;
    }


    protected boolean makeAcceptanceDecision(double currentScore, double newScore){
        if(printDetails()){
            System.out.println("Temperature: " + _temperature);
        }
        if(compareTwoScores(newScore, currentScore)>0)
        {
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
            if(random<acceptanceRatio){
                if(printDetails()){
                    System.out.println("Accept: " + random + " < " + acceptanceRatio);
                    //System.out.println();
                }
                return true;
            }
            else {
                if(printDetails()){
                    System.out.println("Reject: " + random + " > " + acceptanceRatio);
                    //System.out.println();
                }
                return false;
            }
        }
    }

    protected abstract void updateTemperature();

    protected abstract void initializeTemperature();

    protected abstract boolean concludeSearch();
}
