package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAddition;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimberObservableBase;
import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.*;

import java.util.Comparator;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/22/12
 * Time: 1:47 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class AllNeighboursHillClimberBase<G extends Graph<N,E>,N,E,S> extends HillClimberObservableBase<G, S> {
    protected NetworkWholeNeighbourhoodGenerator<G,N,E> _networkGenerator;
    protected Long _maxGenerations;
    protected int _diameterLimit;


    public AllNeighboursHillClimberBase(NetworkWholeNeighbourhoodGenerator<G,N,E> networkGenerator)
    {
        _networkGenerator = networkGenerator;
        
    }


    protected boolean considerNeighborhood(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, final S bestSeenSolutionScore,final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore) {return true;}

    protected boolean considerSameLevelNeighborhood(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, final S bestSeenSolutionScore,final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore, final int diameterLimit) {
        final Ref<Boolean> sawBetterSolution = new Ref<Boolean>(false);
        final Ref<S> currentBestSeenSolutionScore = new Ref<S>(bestSeenSolutionScore);
 
        Func4<G,Integer,E,E,Boolean> rearrangementComputed = getRearrangementComputedListener(solution, getScore, scoreComparator, bestSeenSolutionScore,
                getBetterSolution, newBestScore, sawBetterSolution, currentBestSeenSolutionScore);

       _networkGenerator.generateHorizontalNeighbours(solution, rearrangementComputed, diameterLimit);
       return sawBetterSolution.get();
    }


    protected boolean considerHigherLevelNeighborhood(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, final S bestSeenSolutionScore,final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore) {
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

    protected HillClimbResult<G,S> search(G bestSeenSolution, Func1<G,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore, boolean startSameLevel)
    {
        boolean sameLevel = startSameLevel;
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


    public HillClimbResult<G,S> search(G solution, Func1<G,S> getScore, Comparator<S> scoreComparator, Long maxExaminations, Long maxGenerations, boolean startSameLevel, int diameterLimit)
    {
        _maxExaminations = maxExaminations;
        _maxGenerations = maxGenerations;
        _diameterLimit = diameterLimit;
        S score = getScore.execute(solution);
        this.fireInitialSolutionScoreComputedEvent(solution, score);
        return search(solution, getScore, scoreComparator, score, startSameLevel);
    }

    public HillClimbResult<G,S> search(G solution, Func1<G,S> getScore, Comparator<S> scoreComparator, Long maxExaminations, Long maxGenerations, int diameterLimit)
    {
        _maxExaminations = maxExaminations;
        _maxGenerations = maxGenerations;
        _diameterLimit = diameterLimit;
        S score = getScore.execute(solution);
        this.fireInitialSolutionScoreComputedEvent(solution, score);
        return search(solution, getScore, scoreComparator, score, true);
    }


    public HillClimbResult<G,S> search(G solution, Func1<G,S> getScore, Comparator<S> scoreComparator, Long maxExaminations, Long maxGenerations, boolean startSameLevel)
    {
        return search(solution, getScore, scoreComparator, maxExaminations, maxGenerations, startSameLevel, 0);
    }


    public HillClimbResult<G,S> search(G solution, Func1<G,S> getScore, Comparator<S> scoreComparator, Long maxExaminations, Long maxGenerations)
    {
        return search(solution, getScore, scoreComparator, maxExaminations, maxGenerations, true, 0);
    }



    protected abstract Func4<G,Integer,E,E,Boolean> getRearrangementComputedListener(G solution, Func1<G, S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore,
                                                                               Ref<Func1<G, G>> getBetterSolution, Ref<S> newBestScore, Ref<Boolean> sawBetterSolution,
                                                                               Ref<S> currentBestSeenSolutionScore);

}
