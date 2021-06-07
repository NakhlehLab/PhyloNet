package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.NetworkValidator;
import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func4;

import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/7/12
 * Time: 4:36 PM
 * To change this template use File | Settings | File Templates.
 */
public interface ReticulateEdgeAddition<G,N,E> extends NetworkValidator<N,E>
{
    Map<N,Set<N>> computeRearrangements(G network, boolean validateNetwork, Func4<G,E,E,E,Boolean> rearrangementComputed);

    <S extends Set<N>> void computeRearrangements(G network, boolean validateNetwork, Func4<G,E,E,E,Boolean> rearrangementComputed, Map<N,S> nodeToAncestors);

    G performRearrangement(G network, boolean validateNetwork, E sourceEdge, E destinationEdge);

    <S extends Set<N>> G performRearrangement(G network, boolean validateNetwork, E sourceEdge, E destinationEdge, Map<N,S> nodeToAncestors, Func<S> makeSet);
}
