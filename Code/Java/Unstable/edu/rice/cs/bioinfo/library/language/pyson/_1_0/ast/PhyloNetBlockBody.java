package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 2:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloNetBlockBody implements  Block{

    public final Iterable<PhyloNetCommand> Commands;

    public PhyloNetBlockBody(Iterable<PhyloNetCommand> commands)
    {
        Commands = commands;
    }

    public <R, T, E extends Exception> R execute(BlockAlgo<R, T, E> algo, T input) throws E {
        return algo.forPhylonetBlockBody(this, input);
    }
}
