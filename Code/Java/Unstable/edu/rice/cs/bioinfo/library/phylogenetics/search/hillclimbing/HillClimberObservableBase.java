package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing;
import edu.rice.cs.bioinfo.library.phylogenetics.search.ObservableGenerationalScoringSearcherBase;
import edu.rice.cs.bioinfo.library.programming.*;
import sun.java2d.SunGraphicsEnvironment;

import java.util.Comparator;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 12:19 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class HillClimberObservableBase<T1,S> extends ObservableGenerationalScoringSearcherBase<T1,S> implements HillClimberObservable<T1, S>
{
    private boolean _continueSearch = true;

    private void concludeSearch()
    {
        _continueSearch = false;
    }

    protected  boolean getContinueSearch()
    {
        return _continueSearch;
    }

    private Long _maxExaminations;

    @Override
    public HillClimbResult<T1,S> search(T1 solution, Func1<T1,S> getScore, Comparator<S> scoreComparator)
    {
        _maxExaminations = null;
        S score = getScore.execute(solution);
        this.fireInitialSolutionScoreComputedEvent(solution, score);
        return search(solution, getScore, scoreComparator, score);
    }

    @Override
    public HillClimbResult<T1,S> search(T1 solution, Func1<T1,S> getScore, Comparator<S> scoreComparator, long maxExaminations)
    {
        _maxExaminations = maxExaminations;
        S score = getScore.execute(solution);
       this.fireInitialSolutionScoreComputedEvent(solution, score);
        return search(solution, getScore, scoreComparator, score);
    }


    private HillClimbResult<T1,S> search(T1 bestSeenSolution, Func1<T1,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore)
    {
        for(; _continueSearch; this.incrementGenerationNumber())
        {
            Ref<Func1<T1,T1>>  getBetterNeighbor = new Ref<Func1<T1, T1>>(null);
            Ref<S> newBestScore = new Ref<S>(null);
            boolean sawBetterSolution = considerNeighborhood(bestSeenSolution, getScore, scoreComparator, bestSeenSolutionScore, getBetterNeighbor, newBestScore);
            if(sawBetterSolution)
            {
                bestSeenSolution = getBetterNeighbor.get().execute(bestSeenSolution);
                bestSeenSolutionScore = newBestScore.get();
            }

            if(!sawBetterSolution)
            {
                concludeSearch();
            }
        }

        return new HillClimbResult<T1,S>(bestSeenSolution, bestSeenSolutionScore, getExaminationsCount(), getGenerationNumber());
    }

    protected abstract boolean considerNeighborhood(T1 solution, Func1<T1,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore,
                                                    Ref<Func1<T1,T1>>  getBetterSolution, Ref<S> newBestScore);

    protected boolean considerSolution(T1 solution, Func1<T1,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore, Ref<S> newBestScore)
    {
        incrementExaminations();
        if(_maxExaminations != null && _maxExaminations <= getExaminationsCount())
        {
            concludeSearch();
        }

        S solutionScore =  getScore.execute(solution);
        if(scoreComparator.compare(solutionScore, bestSeenSolutionScore) == 1)
        {
            newBestScore.set(solutionScore);
            this.fireBetterSolutionFoundEvent(solution, solutionScore);
            return true;
        }
        return false;
    }
}
