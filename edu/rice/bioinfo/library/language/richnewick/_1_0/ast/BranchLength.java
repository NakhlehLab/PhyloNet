package edu.rice.bioinfo.library.language.richnewick._1_0.ast;


public interface BranchLength extends  AbstractSyntaxNode{

    public <R,T,E extends Exception> R execute(BranchLengthAlgo<R,T,E> algo, T input) throws E;
}
