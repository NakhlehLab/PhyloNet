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
 *
 * This class is a subclass of SimpleHillClimbing
 * It implements the simulated annealing
 *
 * See "A Maximum Pseudo-likelihood Approach for Phylogenetic Networks", BMC Genomics, 2015.
 */
public abstract class SimulatedAnnealingBase extends SimpleHillClimbing {
    Random _random;
    protected double _temperature = -1;


    /**
     * Constructor / Initializer
     */
    public SimulatedAnnealingBase(Comparator<Double> scoreComparator, NetworkNeighbourhoodGenerator generator, Long seed) {
        super(scoreComparator, generator);
        initializeTemperature();
        if(seed!=null) {
            _random = new Random(seed);
        }else{
            _random = new Random();
        }
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
        updateTemperature();
        return newScore;
    }


    /**
     * This function is to decide whether the proposed species network should be accepted
     */
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
     * This function is to update the temperature
     */
    protected abstract void updateTemperature();


    /**
     * This function is to initialize the temperature
     */
    protected abstract void initializeTemperature();


    /**
     * This function is to decide whether the search should be terminated
     */
    protected abstract boolean concludeSearch();
}
