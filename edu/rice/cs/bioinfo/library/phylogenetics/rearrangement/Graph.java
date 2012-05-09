package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/6/12
 * Time: 2:41 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Graph<N,E> extends GraphReadOnly<N,E>
{
    void addNode(N node);

    void removeNode(N node);

    void addEdge(E edge);

    void removeEdge(E edge);
}
