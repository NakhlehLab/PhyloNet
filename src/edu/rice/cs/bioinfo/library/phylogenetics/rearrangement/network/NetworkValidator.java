package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network;

import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/7/12
 * Time: 5:52 PM
 * To change this template use File | Settings | File Templates.
 */
public interface NetworkValidator<N,E>
{
    void assertValidNetwork(GraphReadOnly<N,E> network);
}
