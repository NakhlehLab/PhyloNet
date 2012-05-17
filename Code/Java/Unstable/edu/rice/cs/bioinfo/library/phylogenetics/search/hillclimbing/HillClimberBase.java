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
public abstract class HillClimberBase<T1,T2> implements HillClimber<T1>
{
    @Override
    public <S> HillClimbResult<T1,S> search(T1 solution, Func1<T1,S> getScore, Comparator<S> scoreComparator)
    {
        S score = getScore.execute(solution);
        return search(solution, getScore, scoreComparator, null, score);
    }

    @Override
    public <S> HillClimbResult<T1,S> search(T1 solution, Func1<T1,S> getScore, Comparator<S> scoreComparator, int maxIterations)
    {
        S score = getScore.execute(solution);
        return search(solution, getScore, scoreComparator, maxIterations, score);
    }


    private <S> HillClimbResult<T1,S> search(T1 bestSeenSolution, Func1<T1,S> getScore, Comparator<S> scoreComparator, Integer maxIterations, S bestSeenSolutionScore)
    {
        Ref<T2> searchState = new Ref<T2>(null);
        for(int i = 1; true; i++)
        {
            Ref<Func1<T1,T1>>  getBetterNeighbor = new Ref<Func1<T1, T1>>(null);
            Ref<S> newBestScore = new Ref<S>(null);
            boolean sawBetterSolution = considerNeighborhood(bestSeenSolution, getScore, scoreComparator, bestSeenSolutionScore, i, getBetterNeighbor, newBestScore, searchState);
            if(sawBetterSolution)
            {
                bestSeenSolution = getBetterNeighbor.get().execute(bestSeenSolution);
                bestSeenSolutionScore = newBestScore.get();
            }

            if((!sawBetterSolution )|| (maxIterations == null ? false : i == maxIterations.intValue()))
            {
                return new HillClimbResult<T1,S>(bestSeenSolution, bestSeenSolutionScore, i);
            }

        }
    }

    protected abstract <S> boolean considerNeighborhood(T1 solution, Func1<T1,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore, int iterationNumber,
                                                        Ref<Func1<T1,T1>>  getBetterSolution, Ref<S> newBestScore, Ref<T2> searchState);

    protected <S> boolean considerSolution(T1 solution, Func1<T1,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore, Ref<S> newBestScore)
    {
        S solutionScore =  getScore.execute(solution);
        if(scoreComparator.compare(solutionScore, bestSeenSolutionScore) == 1)
        {
            newBestScore.set(solutionScore);
            return true;
        }
        return false;
    }
}
