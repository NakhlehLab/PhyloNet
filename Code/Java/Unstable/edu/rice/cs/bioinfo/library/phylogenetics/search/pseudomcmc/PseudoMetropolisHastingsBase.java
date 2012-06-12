package edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc;
import edu.rice.cs.bioinfo.library.phylogenetics.search.ObservableGenerationalScoringSearcherBase;
import edu.rice.cs.bioinfo.library.programming.*;
import org.apache.commons.io.filefilter.FalseFileFilter;

import java.util.Comparator;
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

    private T _bestSeenSolution;

    private S _bestSeenSolutionScore;

    private Func1<T,S> _getScore;

    private Func2<S,S,Double> _divideScore;

    private  boolean _maximize;

    private S _searchParentScore;

    private Random _rand;

      @Override
      public PseudoMetropolisHastingsResult<T, S> search(T solution, Func1<T,S> getScore, Func2<S,S,Double> divideScore, boolean maximize, Random randomSource, long maxExaminations)
      {
          _maxExaminations = maxExaminations;
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


        this.incrementGenerationNumber();
        for(; _continueSearch; this.incrementGenerationNumber())
        {
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

        if(acceptProb >= 1 || _rand.nextDouble() <= acceptProb)
        {
            _bestSeenSolution = solution.DeepCopy();
            _bestSeenSolutionScore = solutionScore;
            this.fireBetterSolutionFoundEvent(solution, solutionScore);
            _searchParentScore = solutionScore;
            return false;
        }
        else
        {
            return true && _continueSearch;
        }

    }
}
