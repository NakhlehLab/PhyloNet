package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 8:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class HillClimbResult<T, S>
{
    public final T LocalOptimum;

    public final S LocalOptimumScore;

    public final int NumIterationsPerformed;

    public HillClimbResult(T localOptimum, S localOptimumScore, int numIterationsPerformed)
    {
        LocalOptimum = localOptimum;
        LocalOptimumScore = localOptimumScore;
        NumIterationsPerformed = numIterationsPerformed;
    }
}
