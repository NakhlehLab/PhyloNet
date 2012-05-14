package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing;

import edu.rice.cs.bioinfo.library.programming.*;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 12:19 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class HillClimberBase<T> implements HillClimber<T>
{
    @Override
    public <S> HillClimbResult<T,S> search(T solution, Func1<T,S> getScore, Comparator<S> scoreComparator)
    {
        S score = getScore.execute(solution);
        return search(solution, getScore, scoreComparator, new Ref<S>(score));
    }


    private <S> HillClimbResult<T,S> search(T bestSeenSolution, Func1<T,S> getScore, Comparator<S> scoreComparator, Ref<S> bestSeenSolutionScore)
    {
        Ref<Func1<T,T>>  getBetterSolution = new Ref<Func1<T, T>>(null);

        boolean sawBetterSolution = considerNeighborhood(bestSeenSolution, getScore, scoreComparator, bestSeenSolutionScore, getBetterSolution);
        if(sawBetterSolution)
        {
            bestSeenSolution = getBetterSolution.get().execute(bestSeenSolution);
            return search(bestSeenSolution, getScore, scoreComparator, bestSeenSolutionScore);
        }
        else
        {
            return new HillClimbResult<T,S>(bestSeenSolution, bestSeenSolutionScore.get());
        }
    }

    protected abstract <S> boolean considerNeighborhood(T solution, Func1<T,S> getScore, Comparator<S> scoreComparator, Ref<S> bestSeenSolutionScore, Ref<Func1<T,T>>  getBetterSolution);

    protected <S> boolean considerSolution(T solution, Func1<T,S> getScore, Comparator<S> scoreComparator, Ref<S> bestSeenSolutionScore)
    {
        S solutionScore =  getScore.execute(solution);
        if(scoreComparator.compare(solutionScore, bestSeenSolutionScore.get()) == 1)
        {
            bestSeenSolutionScore.set(solutionScore);
            return true;
        }
        return false;
    }
}
