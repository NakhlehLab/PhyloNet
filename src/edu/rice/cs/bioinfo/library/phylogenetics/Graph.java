package edu.rice.cs.bioinfo.library.phylogenetics;

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

    void addNodes(N... nodes);

    void removeNode(N node);

    void addEdge(E edge);

    void addEdges(E... edges);

    void removeEdge(E edge);
}
