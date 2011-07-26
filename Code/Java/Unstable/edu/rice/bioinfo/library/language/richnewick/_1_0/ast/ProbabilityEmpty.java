package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

public final class ProbabilityEmpty implements Probability
{
    public static final ProbabilityEmpty Singleton = new ProbabilityEmpty();

    private ProbabilityEmpty()
    {
    }

    public <R, T, E extends Exception> R execute(ProbabilityAlgo<R, T, E> algo, T input) throws E {
        return algo.forProbabilityEmpty(this, input);
    }
}
