package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model;

import jeigen.DenseMatrix;
import java.util.concurrent.ConcurrentHashMap;

/**
 * This abstract class is used to represent nucleotide substitution models.
 * Extending classes only need to implement the necessary methods from the RateModel interface.
 */
public abstract class SubstitutionModel implements RateModel
{

    private ConcurrentHashMap<Double, DenseMatrix> keep = new ConcurrentHashMap<Double, DenseMatrix>();

    /**
     * Gets the probability transition matrix for transitioning between bases after a certain amount of time.
     *
     * @param time The time for transitioning.
     * @return The probability transition matrix.
     */
    final public DenseMatrix getProbabilityMatrix(double time)
    {
        return getRateMatrix().mul(time).mexp();
    }

    /**
     * A cached version of {@link #getProbabilityMatrix(double)}
     * Warning may recompute in case of threads;
     * @param time The time for transitioning.
     * @return The probabilty transition matrix.
     */
    final public DenseMatrix getCachedProbabilityMatrix(double time)
    {
        DenseMatrix result = keep.get(time);
        if (result == null)
        {
            result = getProbabilityMatrix(time);
            keep.put(time,result);
            return result;
        }
        else
            return result;
    }


    public abstract DenseMatrix getProbabilityMatrixIntegrated();

}
