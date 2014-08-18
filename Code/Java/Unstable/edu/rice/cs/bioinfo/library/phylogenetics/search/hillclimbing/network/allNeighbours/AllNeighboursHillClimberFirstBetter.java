package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.*;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Ref;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/26/12
 * Time: 2:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class AllNeighboursHillClimberFirstBetter<G extends Graph<N,E>,N,E,S> extends AllNeighboursHillClimberBase<G,N,E,S> {
    private NetworkNeighbourhoodRandomWalkGenerator<G,N,E> _networkGenerator;
    private long _maxFailure = -1;
    private long _hasTried = 1;

    public AllNeighboursHillClimberFirstBetter(NetworkNeighbourhoodRandomWalkGenerator<G, N, E> generator) {
        _networkGenerator = generator;
    }

    public boolean considerRandomNeighbor(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, final S bestSeenSolutionScore,final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore, final int diameterLimit) {
        final Ref<Boolean> sawBetterSolution = new Ref<Boolean>(false);
        final Ref<S> currentBestSeenSolutionScore = new Ref<S>(bestSeenSolutionScore);

        Func4<G,Integer,E,E,Boolean> rearrangementComputed = getRearrangementComputedListener(solution, getScore, scoreComparator, bestSeenSolutionScore,
                getBetterSolution, newBestScore, sawBetterSolution, currentBestSeenSolutionScore);
        boolean incrementHybrid = getGenerationNumber()==_maxGenerations ? false:true;
        _networkGenerator.computeRandomNeighbour(solution, incrementHybrid, rearrangementComputed, diameterLimit);

        return sawBetterSolution.get();
    }

    public HillClimbResult<G,S> search(G solution, Func1<G,S> getScore, Comparator<S> scoreComparator, Long maxExaminations, int maxGenerations, Long maxFailure, int diameterLimit){
        _maxFailure = maxFailure;
        return super.search(solution, getScore, scoreComparator, maxExaminations, maxGenerations, diameterLimit);
    }

    public HillClimbResult<G,S> search(G solution, Func1<G,S> getScore, Comparator<S> scoreComparator, Long maxExaminations, int maxGenerations, Long maxFailure, int diameterLimit, int hasTried){
        _maxFailure = maxFailure;
        _hasTried = hasTried;
        return super.search(solution, getScore, scoreComparator, maxExaminations, maxGenerations, diameterLimit);
    }



    public HillClimbResult<G,S> search(G bestSeenSolution, Func1<G,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore)
    {
        long tried = _hasTried;
        //System.out.print("Trying ");
        while(_continueSearch)
        {

            if(_maxExaminations != null && _maxExaminations <= getExaminationsCount() || _maxFailure == tried)
            {
                concludeSearch();
            }
            else {
                //System.out.print(" #"+tried);
                Ref<Func1<G, G>> getBetterNeighbor = new Ref<Func1<G, G>>(null);
                Ref<S> newBestScore = new Ref<S>(null);
                boolean sawBetterSolution = considerRandomNeighbor(bestSeenSolution, getScore, scoreComparator, bestSeenSolutionScore, getBetterNeighbor, newBestScore, _diameterLimit);
                if (sawBetterSolution) {
                    tried = 0;
                    //System.out.println("Success: " + bestSeenSolutionScore);
                    bestSeenSolution = getBetterNeighbor.get().execute(bestSeenSolution);
                    bestSeenSolutionScore = newBestScore.get();
                } else {
                    tried++;
                    //System.out.println("Failure #"+tried+": " + bestSeenSolutionScore);
                }
                incrementExaminations();
            }
         }
        //System.out.println();
        return new HillClimbResult<G,S>(bestSeenSolution, bestSeenSolutionScore, getExaminationsCount(), getGenerationNumber());
    }



    @Override
    protected Func4<G,Integer,E,E,Boolean> getRearrangementComputedListener(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, S bestSeenSolutionScore,
                                                                            final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore,
                                                                            final Ref<Boolean> sawBetterSolution, final Ref<S> currentBestSeenSolutionScore)
    {
        return new Func4<G, Integer, E, E, Boolean>() {
            public Boolean execute(G neighbor, final Integer operation, final E edge1, final E edge2)
            {
                boolean relax = operation==1;

                if (considerSolution(neighbor, getScore, scoreComparator, currentBestSeenSolutionScore.get(), newBestScore, relax))
                {
                    sawBetterSolution.set(true);
                    currentBestSeenSolutionScore.set(newBestScore.get());
                    getBetterSolution.set(new Func1<G, G>() {
                        public G execute(G network) {
                            if(operation == 0){
                                incrementGenerationNumber();
                            }
                            else if(operation == 1){
                                decrementGenerationNumber();
                            }
                            return _networkGenerator.performRearrangement(network,operation,edge1,edge2);
                        }
                    });
                    return false;
                }
                return getContinueSearch();
            }
        };
    }

    private boolean considerSolution(G solution, Func1<G,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore, Ref<S> newBestScore, boolean relax)
    {

        S solutionScore =  getScore.execute(solution);

        //System.out.print(solutionScore +" vs. "+ bestSeenSolutionScore);
        if(scoreComparator.compare(solutionScore, bestSeenSolutionScore) == 1)
        {
            //System.out.println(": better");
            newBestScore.set(solutionScore);
            this.fireBetterSolutionFoundEvent(solution, solutionScore);
            return true;
        }
        //System.out.println(": worse");
        return false;
    }
}