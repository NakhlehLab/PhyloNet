package edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;

import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/8/12
 * Time: 11:21 AM
 * To change this template use File | Settings | File Templates.
 */
public interface PseudoMetropolisHastings<T,S>
{
    public PseudoMetropolisHastingsResult<T, S> search(T solution, Func1<T,S> getScore, Func2<S,S,Double> divideScore, boolean maximize, Random randomSource);

    public PseudoMetropolisHastingsResult<T, S> search(T solution, Func1<T,S> getScore, Func2<S,S,Double> divideScore, boolean maximize, Random randomSource, long maxExaminations, long maxGeneration);
}
