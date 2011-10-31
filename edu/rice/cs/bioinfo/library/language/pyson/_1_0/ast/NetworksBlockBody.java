package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/20/11
 * Time: 11:26 AM
 * To change this template use File | Settings | File Templates.
 */
public class NetworksBlockBody extends RNewichAssignmentsBlockBodyBase<NetworkAssignment> {

    public NetworksBlockBody(Iterable<NetworkAssignment> assignments)
    {
        super(assignments);
    }

    public <R, T, E extends Exception> R execute(BlockAlgo<R, T, E> algo, T input) throws E {
        return algo.forNetworksBlock(this, input);
    }
}
