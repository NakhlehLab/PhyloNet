package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

public interface Probability extends AbstractSyntaxNode
{
    public <R,T,E extends Exception> R execute(ProbabilityAlgo<R,T,E> algo, T input) throws E;
}
