package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference;
import java.util.Set;

public abstract class NeighborJoinerTemplate<N,E,G,D extends Comparable<D>> implements NeighborJoiner<N,G>
{
    public G performJoin(Set<N> taxa)
    {
       G resultTree = makeEmptyGraph();

        // complete algo here

       return resultTree;
    }

    /***
     * Returns the distance between two distinct taxa.
     *
     * @param taxon1 a taxon
     * @param taxon2 another taxon
     * @return the distance between teh two taxa.
     */
    protected abstract D getDistance(N taxon1, N taxon2);

    /**
     * Makes a new graph of type G without any member nodes or edges.
     *
     * @return an empty graph.
     */
    protected abstract G makeEmptyGraph();

    /**
     * Creates a new node and adds it to the given graph.
     *
     * @param graph the graph into which to add the returned node
     * @return the node added to the graph.
     */
    protected abstract N makeNewNodeInGraph(G graph);

    /**
     * Creates a edge and adds it to the given graph.
     *
     * @param node1 one of the nodes of the edge
     * @param node2 the other node of the edge
     * @param graph the graph into which to add the new edge
     * @return the created edge
     */
    protected abstract E addEdgeToGraph(N node1, N node2, G graph);

    /**
     * Removes a given edge from the given graph.
     *
     * @param edge the edge to be removed
     * @param graph the graph containing the given edge.
     */
    protected abstract void removeEdgeFromGraph(E edge, G graph);

    /**
     * Returns an instance of type D with the given value.
     *
     * @param value the value to be returned in type D
     * @return the given value represented in type D
     */
    protected abstract D makeD(int value);

    /**
     * Sums two values.
     *
     * @param d1 a value
     * @param d2 a value
     * @return the sum of the given values
     */
    protected abstract D add(D d1, D d2);

    /**
     * Subtracts two values.
     *
     * @param minuend  a value
     * @param subtrahend a value
     * @return
     */
    protected abstract D subtract(D minuend, D subtrahend);

    /**
     * Multiplies two values.
     *
     * @param d1 a value
     * @param d2 a value
     * @return the product of the given distances
     */
    protected abstract D multiply(D d1, D d2);

    /**
     * Divides two values.
     *
     * @param dividend the numerator
     * @param divisor  the denominator
     * @return the division of the given distances.
     */
    protected abstract D divide(D dividend, D divisor);

}
