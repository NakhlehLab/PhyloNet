package edu.rice.cs.bioinfo.programs.phylonet.algos.search.GeneticAlgorithm;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing.HillClimberBase;
import edu.rice.cs.bioinfo.programs.phylonet.algos.search.SearchBase;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;

import java.util.*;

/**
 * Created by yunyu on 11/6/14.
 *
 * This class is a subclass of searchBase
 * It implements the genetic algorithm for search
 */
public abstract class GeneticAlgorithmBase extends SearchBase{
    private int _maxGenerations;
    private int _generationCounter = 1;


    /**
     * Constructor / Initializer
     */
    public GeneticAlgorithmBase(Comparator<Double> scoreComparator){
        super(scoreComparator);

    }


    /**
     * This function is to reset the number of generations
     */
    protected void clearGenerationCounter(){
        _generationCounter = 1;
    }



    /**
     * This function is to search the network space for one round
     *
     * @param startNetworks             the starting generations
     * @param getScore                  the function that computes the score of a candidate species network
     * @param numOptimums               the number of optimal species networks to return
     * @param maxGenerations            the maximal number of generations for the search; served as stopping criterion
     * @param scoreEachTopologyOnce     indicate whether each network topology only need to be evaluated once
     * @param resultList                the resulting species networks along with their scores
     */
    public void search(List<Network> startNetworks, Func1<Network,Double> getScore, int numOptimums, int maxGenerations, boolean scoreEachTopologyOnce, LinkedList<Tuple<Network,Double>> resultList){
        _maxGenerations = maxGenerations;
        _optimalNetworks = resultList;
        _numOptimums = numOptimums;
        _scoreEachTopologyOnce = scoreEachTopologyOnce;
        if(getCashedNetworkListSize()!=0 && _optimalNetworks.size()==0){
            for(Tuple<Network,Double> tried: getAllCashedNetworks()){
                updateOptimalNetworks(tried.Item1, tried.Item2);
            }
        }
        LinkedList<MutableTuple<Network,Double>> startingGeneration = new LinkedList<>();
        for(Network network: startNetworks){
            startingGeneration.add(new MutableTuple<Network, Double>(network, 0.0));
        }
        computeScoresForOneGeneration(startingGeneration, getScore);
        search(startingGeneration, getScore);
    }



    /**
     * This is the main function to search the network space for one round
     *
     * @param currentGeneration         the current generation along with their scores
     * @param getScore                  the function that computes the score of a candidate species network     the resulting species networks along with their scores
     */
    protected void search(LinkedList<MutableTuple<Network,Double>> currentGeneration, Func1<Network,Double> getScore)
    {
        while(!concludeSearch())
        {
            if(printDetails()){
                System.out.println("Generation #" + _generationCounter+":");
                for(MutableTuple<Network,Double> individual: currentGeneration){
                    System.out.println(individual.Item2+":"+individual.Item1.toString());
                }
                System.out.println();
            }

            LinkedList<MutableTuple<Network,Double>> nextGeneration = generateNextGeneration(currentGeneration);
            computeScoresForOneGeneration(nextGeneration, getScore);
            updateResults(nextGeneration);
            currentGeneration = nextGeneration;
            _generationCounter++;
        }
    }


    /**
     * This function is to update the caches and results
     *
     * @param currentGeneration     the current generation along with their scores
     */
    private void updateResults(LinkedList<MutableTuple<Network,Double>> currentGeneration){
        for(MutableTuple<Network,Double> individual: currentGeneration){
            super.updateCashedResults(individual.Item1, individual.Item2);
            super.updateOptimalNetworks(individual.Item1, individual.Item2);
        }
    }


    /**
     * This function is to generate the next generation based on the current generation
     *
     * @param currentGeneration     the current generation along with their scores
     */
    protected abstract LinkedList<MutableTuple<Network,Double>> generateNextGeneration(LinkedList<MutableTuple<Network,Double>> currentGeneration);



    /**
     * This function is to compute the scores for each of the species network in the current generation
     *
     * @param currentGeneration     the current generation along with their scores
     */
    protected void computeScoresForOneGeneration(LinkedList<MutableTuple<Network,Double>> currentGeneration, Func1<Network,Double> getScore){
        for(MutableTuple<Network,Double> individual: currentGeneration){
            individual.Item2 = getScore.execute(individual.Item1);
        }
    }


    /**
     * This function is to decide whether the search should be terminated
     */
    protected boolean concludeSearch(){
        return _generationCounter >= _maxGenerations;
    }

}
