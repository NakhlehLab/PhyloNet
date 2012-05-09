package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.*;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.*;
import edu.rice.cs.bioinfo.library.programming.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/7/12
 * Time: 4:36 PM
 * To change this template use File | Settings | File Templates.
 */
public interface ReticulateEdgeAddition<G,N,E> extends NetworkValidator<N,E>
{
    void computeRearrangementsWithValidation(G network, Proc4<G,E,E,E> rearrangementComputed);

    void computeRearrangementsWithoutValidation(G network, Proc4<G,E,E,E> rearrangementComputed);
}
