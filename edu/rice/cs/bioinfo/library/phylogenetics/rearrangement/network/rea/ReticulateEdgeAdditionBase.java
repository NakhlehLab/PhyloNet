package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.NetworkValidatorBase;
import edu.rice.cs.bioinfo.library.programming.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/7/12
 * Time: 5:47 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReticulateEdgeAdditionBase<G extends GraphReadOnly<N,E>,N,E> extends NetworkValidatorBase<N,E> implements ReticulateEdgeAddition<G,N,E>
{
    public void computeRearrangementsWithValidation(G network, Proc4<G,E,E,E> rearrangementComputed)
    {
        this.assertValidNetwork(network);
        computeRearrangementsWithoutValidation(network, rearrangementComputed);
    }
}
