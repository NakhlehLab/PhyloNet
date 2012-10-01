package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network;

import edu.rice.cs.bioinfo.library.phylogenetics.GraphValidator;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/7/12
 * Time: 5:53 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkValidatorBase<N,E> implements NetworkValidator<N,E>
{
    public void assertValidNetwork(GraphReadOnly<N, E> network)
    {
        new GraphValidator().assertValidGraph(network);
    }
}
