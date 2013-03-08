package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkWholeNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Ref;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/22/12
 * Time: 2:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class AllNeighboursHillClimberSteepestAscent <G extends Graph<N,E>,N,E,S> extends AllNeighboursHillClimberBase<G,N,E,S> {
    private NetworkWholeNeighbourhoodGenerator<G,N,E> _networkGenerator;

    public AllNeighboursHillClimberSteepestAscent(NetworkWholeNeighbourhoodGenerator<G, N, E> generator) {
        _networkGenerator = generator;
    }


    public boolean considerSameLevelNeighborhood(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, final S bestSeenSolutionScore,final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore, final int diameterLimit) {
        final Ref<Boolean> sawBetterSolution = new Ref<Boolean>(false);
        final Ref<S> currentBestSeenSolutionScore = new Ref<S>(bestSeenSolutionScore);

        Func4<G,Integer,E,E,Boolean> rearrangementComputed = getRearrangementComputedListener(solution, getScore, scoreComparator, bestSeenSolutionScore,
                getBetterSolution, newBestScore, sawBetterSolution, currentBestSeenSolutionScore);

        _networkGenerator.generateHorizontalNeighbours(solution, rearrangementComputed, diameterLimit);
        return sawBetterSolution.get();
    }


    public boolean considerHigherLevelNeighborhood(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, final S bestSeenSolutionScore,final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore) {
        final Ref<Boolean> sawBetterSolution = new Ref<Boolean>(false);
        final Ref<S> currentBestSeenSolutionScore = new Ref<S>(bestSeenSolutionScore);

        Func4<G,Integer,E,E,Boolean> rearrangementComputed = getRearrangementComputedListener(solution, getScore, scoreComparator, bestSeenSolutionScore,
                getBetterSolution, newBestScore, sawBetterSolution, currentBestSeenSolutionScore);

        _networkGenerator.generateVerticalNeighbours(solution, rearrangementComputed);
        return sawBetterSolution.get();
    }

/*
    protected HillClimbResult<G,S> search(G bestSeenSolution, Func1<G,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore)
    {
        boolean sameLevel = true;
        while(_continueSearch)
        {
            //System.out.println("start searching same level : " + sameLevel);
            Ref<Func1<G,G>>  getBetterNeighbor = new Ref<Func1<G, G>>(null);
            Ref<S> newBestScore = new Ref<S>(null);
            boolean sawBetterSolution;
            if(sameLevel){
                sawBetterSolution =considerSameLevelNeighborhood(bestSeenSolution, getScore, scoreComparator, bestSeenSolutionScore, getBetterNeighbor, newBestScore);
            }
            else{
                sawBetterSolution =considerHigherLevelNeighborhood(bestSeenSolution, getScore, scoreComparator, bestSeenSolutionScore, getBetterNeighbor, newBestScore);
                this.incrementGenerationNumber();
            }
            if(sawBetterSolution)
            {
                bestSeenSolution = getBetterNeighbor.get().execute(bestSeenSolution);
                bestSeenSolutionScore = newBestScore.get();
                if(!sameLevel){
                    sameLevel = true;
                }
            }
            if(!sawBetterSolution )
            {
                if(!sameLevel){
                    concludeSearch();
                }else{
                    if(getGenerationNumber()==_maxGenerations){
                        concludeSearch();
                    }
                    else{
                        sameLevel = false;
                    }
                }
            }
            //System.out.println(sawBetterSolution + ": " + bestSeenSolutionScore);
        }

        return new HillClimbResult<G,S>(bestSeenSolution, bestSeenSolutionScore, getExaminationsCount(), getGenerationNumber());
    }
    */

    public HillClimbResult<G,S> search(G bestSeenSolution, Func1<G,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore)
    {
        boolean sameLevel = true;
        while(_continueSearch)
        {

            Ref<Func1<G,G>>  getBetterNeighbor = new Ref<Func1<G, G>>(null);
            Ref<S> newBestScore = new Ref<S>(null);
            boolean sawBetterSolution;
            if(sameLevel){
                //System.out.println("start searching same level : " + sameLevel);
                sawBetterSolution =considerSameLevelNeighborhood(bestSeenSolution, getScore, scoreComparator, bestSeenSolutionScore, getBetterNeighbor, newBestScore, _diameterLimit);
            }
            else{
                //System.out.println("start searching higher level : " + sameLevel);
                sawBetterSolution =considerHigherLevelNeighborhood(bestSeenSolution, getScore, scoreComparator, bestSeenSolutionScore, getBetterNeighbor, newBestScore);
                this.incrementGenerationNumber();
            }
            if(sawBetterSolution)
            {
                bestSeenSolution = getBetterNeighbor.get().execute(bestSeenSolution);
                bestSeenSolutionScore = newBestScore.get();
                if(!sameLevel){
                    sameLevel = true;
                }
            }
            if(!sawBetterSolution )
            {
                if(!sameLevel){
                    concludeSearch();
                }else{
                    if(getGenerationNumber()==_maxGenerations){
                        concludeSearch();
                    }
                    else{
                        sameLevel = false;
                    }
                }
            }
            //System.out.println(sawBetterSolution + ": " + bestSeenSolutionScore);
        }
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
                if (AllNeighboursHillClimberSteepestAscent.this.considerSolution(neighbor, getScore, scoreComparator, currentBestSeenSolutionScore.get(), newBestScore))
                {
                    sawBetterSolution.set(true);
                    currentBestSeenSolutionScore.set(newBestScore.get());
                    getBetterSolution.set(new Func1<G, G>() {
                        public G execute(G network) {
                            return _networkGenerator.performRearrangement(network,operation,edge1,edge2);
                        }
                    });
                }
                return getContinueSearch();
            }
        };
    }
}
