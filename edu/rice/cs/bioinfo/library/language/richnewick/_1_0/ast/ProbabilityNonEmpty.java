package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

public class ProbabilityNonEmpty implements Probability
{
    public final Text ProbabilityValue;

    public ProbabilityNonEmpty(Text probabilityValue)
    {
        ProbabilityValue = probabilityValue;
    }


    public <R, T, E extends Exception> R execute(ProbabilityAlgo<R, T, E> algo, T input) throws E {
        return algo.forProbabilityNonEmpty(this, input);
    }
}
