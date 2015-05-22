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
 */
public abstract class GeneticAlgorithmBase extends SearchBase{
    private int _maxGenerations;
    private int _generationCounter = 1;

    public GeneticAlgorithmBase(Comparator<Double> scoreComparator){
        super(scoreComparator);

    }

    protected void clearGenerationCounter(){
        _generationCounter = 1;
    }


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




    protected void search(LinkedList<MutableTuple<Network,Double>> currentGeneration, Func1<Network,Double> getScore)
    {
        while(!concludeSearch())
        {
            //System.out.println("Generation #" + _generationCounter+":");
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

    private void updateResults(LinkedList<MutableTuple<Network,Double>> currentGeneration){
        for(MutableTuple<Network,Double> individual: currentGeneration){
            super.updateCashedResults(individual.Item1, individual.Item2);
            super.updateOptimalNetworks(individual.Item1, individual.Item2);
        }
    }

    protected abstract LinkedList<MutableTuple<Network,Double>> generateNextGeneration(LinkedList<MutableTuple<Network,Double>> currentGenerations);

    /*
    protected abstract LinkedList<MutableTuple<Network,Double>> selectIndividualsForNextGenerations(LinkedList<MutableTuple<Network,Double>> currentGenerations);

    protected abstract void doMutations(List<MutableTuple<Network,Double>> currentGenerations);

    protected abstract void doRecombinations(List<MutableTuple<Network,Double>> currentGenerations);
    */

    protected void computeScoresForOneGeneration(LinkedList<MutableTuple<Network,Double>> currentGeneration, Func1<Network,Double> getScore){
        for(MutableTuple<Network,Double> individual: currentGeneration){
            individual.Item2 = getScore.execute(individual.Item1);
        }
        /*
        //Sorting all individuals
        for(int i=0; i<currentGeneration.size(); i++){
            MutableTuple<Network,Double> individual1 = currentGeneration.get(i);
            for(int j=0;j<i;j++){
                MutableTuple<Network,Double> individual2 = currentGeneration.get(j);
                if(compareTwoScores(individual1.Item2, individual2.Item2)>0){
                    currentGeneration.remove(individual1);
                    currentGeneration.add(j,individual1);
                    break;
                }
            }
        }
        */

    }

    protected boolean concludeSearch(){
        return _generationCounter>=_maxGenerations;
    }

}
