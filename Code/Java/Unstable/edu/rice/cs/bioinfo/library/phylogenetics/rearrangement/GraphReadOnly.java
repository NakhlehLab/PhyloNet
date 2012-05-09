package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement;

import edu.rice.cs.bioinfo.library.programming.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/6/12
 * Time: 2:25 PM
 * To change this template use File | Settings | File Templates.
 */
public interface GraphReadOnly<N,E>
{
    Iterable<N> getNodes();

    Iterable<E> getEdges();

    Tuple<N,N> getNodesOfEdge(E edge);

    Iterable<E> getIncidentEdges(N node);

    boolean isRooted();

    boolean isDestinationNode(N node, E edge);
}
