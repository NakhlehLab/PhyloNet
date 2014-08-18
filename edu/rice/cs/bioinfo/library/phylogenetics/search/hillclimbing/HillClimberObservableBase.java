package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing;

import edu.rice.cs.bioinfo.library.phylogenetics.search.ObservableGenerationalScoringSearcherBase;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Ref;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 12:19 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class HillClimberObservableBase<T1,S> extends ObservableGenerationalScoringSearcherBase<T1,S> implements HillClimberObservable<T1, S>
{
    protected boolean _continueSearch = true;

    protected void concludeSearch()
    {
        _continueSearch = false;
    }

    protected  boolean getContinueSearch()
    {
        return _continueSearch;
    }

    protected Long _maxExaminations;
    protected int _maxGenerations;


    @Override
    public HillClimbResult<T1,S> search(T1 solution, Func1<T1,S> getScore, Comparator<S> scoreComparator)
    {
        _maxExaminations = null;
        S score = getScore.execute(solution);
        this.fireInitialSolutionScoreComputedEvent(solution, score);
        return search(solution, getScore, scoreComparator, score);
    }

    //@Override
    public HillClimbResult<T1,S> search(T1 solution, Func1<T1,S> getScore, Comparator<S> scoreComparator, Long maxExaminations, int maxGenerations)
    {
        _maxExaminations = maxExaminations;
        _maxGenerations = maxGenerations;
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


    protected HillClimbResult<T1,S> search(T1 bestSeenSolution, Func1<T1,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore)
    {
        while(_continueSearch)
        {
            this.incrementGenerationNumber();
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
            //System.out.println(sawBetterSolution + ": " + bestSeenSolutionScore);
            if(_maxGenerations!=-1 && _maxGenerations==getGenerationNumber())break;
            //System.out.println(sawBetterSolution + ": " + bestSeenSolutionScore);
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
