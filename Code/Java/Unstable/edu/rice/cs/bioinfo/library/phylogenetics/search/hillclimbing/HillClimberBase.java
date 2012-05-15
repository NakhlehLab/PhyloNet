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
    public <S> HillClimbResult<T,S> search(T solution, Func1<T,S> getScore, Comparator<S> scoreComparator, int maxIterations)
    {
        S score = getScore.execute(solution);
        return search(solution, getScore, scoreComparator, maxIterations, score);
    }


    private <S> HillClimbResult<T,S> search(T bestSeenSolution, Func1<T,S> getScore, Comparator<S> scoreComparator, int maxIterations, S bestSeenSolutionScore)
    {
        for(int i = 1; true; i++)
        {
            Ref<Func1<T,T>>  getBetterNeighbor = new Ref<Func1<T, T>>(null);
            Ref<S> newBestScore = new Ref<S>(null);
            boolean sawBetterSolution = considerNeighborhood(bestSeenSolution, getScore, scoreComparator, bestSeenSolutionScore, getBetterNeighbor, newBestScore);
            if(sawBetterSolution)
            {
                bestSeenSolution = getBetterNeighbor.get().execute(bestSeenSolution);
                bestSeenSolutionScore = newBestScore.get();
            }

            if(!sawBetterSolution || i == maxIterations)
            {
                return new HillClimbResult<T,S>(bestSeenSolution, bestSeenSolutionScore);
            }

        }
    }

    protected abstract <S> boolean considerNeighborhood(T solution, Func1<T,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore,
                                                        Ref<Func1<T,T>>  getBetterSolution, Ref<S> newBestScore);

    protected <S> boolean considerSolution(T solution, Func1<T,S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore, Ref<S> newBestScore)
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
