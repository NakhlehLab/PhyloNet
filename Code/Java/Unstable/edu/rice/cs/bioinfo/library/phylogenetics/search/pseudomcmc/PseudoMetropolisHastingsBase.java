package edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc;

import edu.rice.cs.bioinfo.library.phylogenetics.search.ObservableGenerationalScoringSearcherBase;
import edu.rice.cs.bioinfo.library.programming.DeepCopyable;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;

import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/8/12
 * Time: 11:26 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class PseudoMetropolisHastingsBase<T extends DeepCopyable<T>,S> extends ObservableGenerationalScoringSearcherBase<T,S> implements PseudoMetropolisHastings<T,S>
{
    protected boolean _continueSearch = true;

    private void concludeSearch()
    {
        _continueSearch = false;
    }

    protected  boolean getContinueSearch()
    {
        return _continueSearch;
    }

    protected Long _maxExaminations = null;

    protected Long _maxGenerations = null;

    protected T _bestSeenSolution;

    protected S _bestSeenSolutionScore;

    protected Func1<T,S> _getScore;

    protected Func2<S,S,Double> _divideScore;

    protected  boolean _maximize;

    protected S _searchParentScore;

    protected Random _rand;

      @Override
      public PseudoMetropolisHastingsResult<T, S> search(T solution, Func1<T,S> getScore, Func2<S,S,Double> divideScore, boolean maximize, Random randomSource, long maxExaminations, long maxGeneration)
      {
          _maxExaminations = maxExaminations;
          _maxGenerations = maxGeneration;
          return search(solution, getScore, divideScore, maximize, randomSource);
      }


    @Override
    public PseudoMetropolisHastingsResult<T, S> search(T solution, Func1<T,S> getScore, Func2<S,S,Double> divideScore, boolean maximize, Random randomSource)
    {
        _getScore = getScore;
        _maximize = maximize;
        _divideScore = divideScore;
        _rand = randomSource;
        _bestSeenSolution = solution.DeepCopy();
        _bestSeenSolutionScore = getScore.execute(solution);
        this.incrementExaminations();
        this.fireInitialSolutionScoreComputedEvent(_bestSeenSolution, _bestSeenSolutionScore);
        _searchParentScore = _bestSeenSolutionScore;


        while(_continueSearch && (_maxGenerations == null ? true : getGenerationNumber() + 1 <= _maxGenerations))
        {
            this.incrementGenerationNumber();
            solution = considerNeighborhood(solution);
        }

        return new PseudoMetropolisHastingsResult<T, S>(_bestSeenSolution, _bestSeenSolutionScore, getExaminationsCount(), getGenerationNumber());
    }


    protected abstract T considerNeighborhood(T currentSolution);

    protected boolean considerSolution(T solution)
    {
        incrementExaminations();
        if(_maxExaminations != null && _maxExaminations <= getExaminationsCount())
        {
            concludeSearch();
        }

        S solutionScore =  _getScore.execute(solution);


        double acceptProb = _maximize ? _divideScore.execute(solutionScore, _searchParentScore):
                                        _divideScore.execute(_searchParentScore, solutionScore);

        boolean solutionIsNewBestSeen = (_maximize ? _divideScore.execute(solutionScore, _bestSeenSolutionScore):
                _divideScore.execute(_bestSeenSolutionScore, solutionScore)) > 1.0;

        if(solutionIsNewBestSeen)
        {
            _bestSeenSolution = solution.DeepCopy();
            _bestSeenSolutionScore = solutionScore;
            this.fireBetterSolutionFoundEvent(solution, solutionScore);
        }

        if(_rand.nextDouble() <= acceptProb)
        {

            _searchParentScore = solutionScore;
            return false;
        }
        else
        {
            return _continueSearch;
        }

    }
}
