package edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc;

import edu.rice.cs.bioinfo.library.programming.DeepCopyable;
import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;

import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/13/12
 * Time: 2:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class PMHGenerationLimitRestartSearcher<T extends DeepCopyable<T>,S> implements PseudoMetropolisHastings<T,S>
{
    private Func<PseudoMetropolisHastings<T,S>> _makeSearcher;

    public PMHGenerationLimitRestartSearcher(Func<PseudoMetropolisHastings<T,S>> makeSearcher)
    {
        _makeSearcher = makeSearcher;
    }

    public PseudoMetropolisHastingsResult<T, S> search(T solution, Func1<T,S> getScore, Func2<S,S,Double> divideScore, boolean maximize, Random randomSource)
    {
        return _makeSearcher.execute().search(solution, getScore, divideScore, maximize, randomSource);
    }

    public PseudoMetropolisHastingsResult<T, S> search(T solution, Func1<T,S> getScore, Func2<S,S,Double> divideScore, boolean maximize, Random randomSource, long maxExaminations, long maxGenerations)
    {
        PseudoMetropolisHastingsResult<T, S> result = _makeSearcher.execute().search(solution.DeepCopy(), getScore, divideScore, maximize, randomSource, maxExaminations, maxGenerations);
        PseudoMetropolisHastingsResult<T, S> bestSeenResult = result;
        long totalExaminations = result.ExaminationsCount;
        long largestGeneration = result.GenerationCount;

        while(totalExaminations < maxExaminations)
        {
            result = _makeSearcher.execute().search(solution.DeepCopy(), getScore, divideScore, maximize, randomSource, maxExaminations - totalExaminations, maxGenerations);
            totalExaminations += result.ExaminationsCount;
            largestGeneration = result.GenerationCount > largestGeneration ? result.GenerationCount : largestGeneration;

            Double ratio = divideScore.execute(bestSeenResult.BestExaminedScore, result.BestExaminedScore);

            if(maximize ? ratio < 1 : ratio > 1)
            {
                bestSeenResult = result;
            }
        }

        return new PseudoMetropolisHastingsResult<T, S>(bestSeenResult.BestExaminedNetwork, bestSeenResult.BestExaminedScore, totalExaminations, largestGeneration);


    }
}
