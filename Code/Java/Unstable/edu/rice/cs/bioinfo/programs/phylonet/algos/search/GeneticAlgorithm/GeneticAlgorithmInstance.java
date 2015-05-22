package edu.rice.cs.bioinfo.programs.phylonet.algos.search.GeneticAlgorithm;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by yunyu on 11/6/14.
 */
public class GeneticAlgorithmInstance extends GeneticAlgorithmBase{
    private double _topologyMutationProb;
    private double _recombinationProb;
    private double _parameterMutationProportion;
    private int _k;
    private NetworkRandomParameterNeighbourGenerator _parametersMutationMaker;
    private NetworkRandomTopologyNeighbourGenerator _topologyMutationMaker;
    private TwoNetworkRandomPDGGenerator _recombinationMaker;
    private double[] _leaveOffspringProb;
    private Random _random;
    private int _generationSize;

    /*public GeneticAlgorithmInstance(Comparator<Double> scoreComparator, double topologyMutationProb, double recombinationProb, double parameterMutationProportion){
        super(scoreComparator);
        _topologyMutationProb = topologyMutationProb;
        _recombinationProb = recombinationProb;
        _parameterMutationProportion = parameterMutationProportion;
        //_kOverN = kOverN;
    }
    */

    public GeneticAlgorithmInstance(int generationSize, Comparator<Double> scoreComparator, NetworkRandomParameterNeighbourGenerator parametersMutator,double parameterMutationProportion, NetworkRandomTopologyNeighbourGenerator topologyMutator, double topologyMutationProb, TwoNetworkRandomPDGGenerator recombiner, double recombinationProb, int k, Long seed){
        super(scoreComparator);
        _generationSize = generationSize;
        _topologyMutationProb = topologyMutationProb;
        _recombinationProb = recombinationProb;
        _parameterMutationProportion = parameterMutationProportion;
        parametersMutator.setParametersToChange(new Func1<Integer, Integer>() {
            @Override
            public Integer execute(Integer input) {
                //return 1;
                return Math.max(1,(int)(input*_parameterMutationProportion));
            }
        });
        //_kOverN = kOverN;
        _parametersMutationMaker = parametersMutator;
        _topologyMutationMaker = topologyMutator;
        _recombinationMaker = recombiner;
        _k = k;
        if(seed==null){
            _random = new Random();
        }
        else{
            _random = new Random(seed);
        }
        double p = 2.0/(generationSize*(generationSize+1));
        _leaveOffspringProb = new double[generationSize];
        int shares = generationSize;
        for(int i=0; i<generationSize; i++){
            if(i==generationSize-1){
                _leaveOffspringProb[i] = 1;
            }
            _leaveOffspringProb[i] = p*(shares--);
            if(i!=0){
                _leaveOffspringProb[i] += _leaveOffspringProb[i-1];
            }
        }

    }

    public GeneticAlgorithmInstance(int generationSize, Comparator<Double> scoreComparator, NetworkRandomParameterNeighbourGenerator parametersMutator, double parameterMutationProportion, NetworkRandomTopologyNeighbourGenerator topologyMutator, double topologyMutationProb, TwoNetworkRandomPDGGenerator recombiner, double recombinationProb, int k){
        this(generationSize,scoreComparator, parametersMutator, parameterMutationProportion, topologyMutator, topologyMutationProb, recombiner, recombinationProb, k, null);
    }

    private LinkedList<MutableTuple<Network,Double>> selectIndividualsForNextGenerationsUsingRanking(LinkedList<MutableTuple<Network,Double>> currentGenerations){
        LinkedList<MutableTuple<Network,Double>> nextGeneration = new LinkedList<>();
        MutableTuple<Network,Double> bestNetwork = currentGenerations.getFirst();
        int n = currentGenerations.size();
        for(int i=0; i<_k; i++){
            nextGeneration.add(new MutableTuple(bestNetwork.Item1.clone(),0));
        }
        for(int i=_k; i<n; i++){
            nextGeneration.add(new MutableTuple(currentGenerations.get(getRandomParentID()).Item1.clone(),0));
        }
        return nextGeneration;
    }


    private LinkedList<MutableTuple<Network,Double>> selectIndividualsForNextGenerationsUsingFitness(LinkedList<MutableTuple<Network,Double>> currentGenerations){
        LinkedList<MutableTuple<Network,Double>> nextGeneration = new LinkedList<>();
        MutableTuple<Network,Double> bestNetwork = currentGenerations.getFirst();
        int n = currentGenerations.size();
        for(int i=0; i<_k; i++){
            nextGeneration.add(new MutableTuple(bestNetwork.Item1.clone(),0));
        }
        int index = 0;
        double total = 0;
        double worst = currentGenerations.getLast().Item2;
        for(MutableTuple<Network,Double> individual: currentGenerations){
            _leaveOffspringProb[index] = individual.Item2-worst;
            total += _leaveOffspringProb[index];
            index++;
        }
        for(int i=0; i<_leaveOffspringProb.length; i++){
            if(i == _leaveOffspringProb.length-1){
                _leaveOffspringProb[i] = 1;
                break;
            }
            _leaveOffspringProb[i] = _leaveOffspringProb[i]/total;
            if(i!=0){
                _leaveOffspringProb[i] +=  _leaveOffspringProb[i-1];
            }

        }
        for(int i=_k; i<n; i++){
            nextGeneration.add(new MutableTuple(currentGenerations.get(getRandomParentID()).Item1.clone(),0));
        }
        return nextGeneration;
    }

    private int getRandomParentID(){
        double random = _random.nextDouble();
        int parentID = 0;
        for(; parentID<_leaveOffspringProb.length; parentID++){
            if(random<_leaveOffspringProb[parentID]){
                break;
            }
        }
        return parentID;
    }


    private int getRandomParentID(int excludeID){
        int parentID;
        do{
            parentID = getRandomParentID();
        }while(parentID!=excludeID);
        return parentID;
    }


    protected LinkedList<MutableTuple<Network,Double>> generateNextGeneration(LinkedList<MutableTuple<Network,Double>> currentGenerations){
        LinkedList<MutableTuple<Network,Double>> nextGenerations = selectIndividualsForNextGenerationsUsingRanking(currentGenerations);
        doMutations(nextGenerations);
        doRecombinations(currentGenerations, nextGenerations);
        return nextGenerations;
    }


    private void doMutations(List<MutableTuple<Network,Double>> currentGenerations) {
        boolean first = true;
        for(MutableTuple<Network,Double> individual: currentGenerations){
            if(first){
                first = false;
                continue;
            }
            //System.out.print("Parameter change: ");
            _parametersMutationMaker.mutateNetwork(individual.Item1);
            //System.out.println("Done");

            if(_random.nextDouble()<_topologyMutationProb){
                //System.out.print("Topology change: ");
                _topologyMutationMaker.mutateNetwork(individual.Item1);
                //System.out.println("Done");
            }
            else{
                //_parametersMutationMaker.mutateNetwork(individual.Item1);
            }
        }
    }


    private void doRecombinations(List<MutableTuple<Network,Double>> parentGenerations, List<MutableTuple<Network,Double>> currentGenerations) {
        int index = 0;
        for(MutableTuple<Network,Double> individual: currentGenerations) {
            if (index != 0) {
                if (_random.nextDouble() < _recombinationProb) {
                    //System.out.print("recombination: ");
                    Network newOffspring = parentGenerations.get(getRandomParentID(index)).Item1.clone();
                    Networks.autoLabelNodes(individual.Item1);
                    Networks.autoLabelNodes(newOffspring);
                    //System.out.println("Parent #1: " + individual.Item1.toString());
                    //System.out.println("Parent #2: " + newOffspring.toString());
                    boolean success = _recombinationMaker.mutateNetwork(individual.Item1, newOffspring);
                    if(success){
                        individual.Item1 = newOffspring;
                        //System.out.println("Offspring: " + newOffspring.toString());
                    }
                    else{
                        //System.out.println("Cannot do recombination!\n");
                    }
                    //System.out.println("Done");
                }
            }
            index++;
        }
    }




    protected void computeScoresForOneGeneration(LinkedList<MutableTuple<Network,Double>> currentGeneration, Func1<Network,Double> getScore){
        super.computeScoresForOneGeneration(currentGeneration, getScore);
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
    }


    public void search(Network startNetwork, Func1<Network,Double> getScore, int numOptimums, int numRuns, int maxGenerations, boolean scoreEachTopologyOnce, LinkedList<Tuple<Network,Double>> resultList){
        for(int i=0; i<numRuns; i++) {
            List<Network> initialGeneration = initializeInitialGeneration(startNetwork.clone());
            clearGenerationCounter();
            super.search(initialGeneration, getScore, numOptimums, maxGenerations, scoreEachTopologyOnce, resultList);
        }
    }

    private List<Network> initializeInitialGeneration(Network startNetwork){

        List<Network> initialGeneration = new ArrayList<>();
        initialGeneration.add(startNetwork);
        double[] threshold = {_generationSize*0.4, _generationSize*0.3, _generationSize*0.2,  _generationSize*0.1};
        int move = 1;
        int count = 1;
        boolean stop = false;
        do{
            for (int i = 0; i < threshold[i]; i++) {
                if(count==_generationSize){
                    stop = true;
                    break;
                }
                Network copy = startNetwork.clone();
                for(int j=0; j<move; j++) {
                    if (_random.nextDouble() < 0.5) {
                        _topologyMutationMaker.mutateNetwork(copy);
                    } else {
                        _parametersMutationMaker.mutateNetwork(copy);
                    }
                }
                initialGeneration.add(copy);
                count++;
            }
            move++;
        }while(!stop);
        return initialGeneration;
    }
}
