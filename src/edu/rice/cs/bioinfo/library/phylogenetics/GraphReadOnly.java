package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Tuple;

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

    E getEdge(N source, N destination);

    boolean isRooted();

    boolean containsEdge(N source, N destination);
}
