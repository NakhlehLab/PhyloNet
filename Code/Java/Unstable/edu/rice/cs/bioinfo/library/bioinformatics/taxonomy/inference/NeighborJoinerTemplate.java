package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference;
import com.sun.xml.internal.xsom.impl.scd.Axis;

import java.util.HashSet;
import java.util.Set;

public abstract class NeighborJoinerTemplate<N,E,G,D extends Comparable<D>> implements NeighborJoiner<N,G>
{
    public D MAX;

    public G performJoin(Set<N> taxa)
    {
        G resultTree = makeEmptyGraph();
        Set<N> Node = new HashSet<N>();
        Node.addAll(taxa);

        while(Node.size() > 1 ) {
            /* find the minimum distance and the corresponding nodes */
            D temp = MAX;
            N node1 = makeNewNodeInGraph(resultTree);
            N node2 = makeNewNodeInGraph(resultTree);
            for(N a: Node){
                for(N b: Node){
                    if(!a.equals(b)){
                        if(temp.compareTo(getDistance(a, b)) > 0){
                            temp = getDistance(a, b);
                            node1 = a;
                            node2 = b;
                        }
                    }
                }
            }
            /* delete node1 and node2 from Set */
            Node.remove(node1);
            Node.remove(node2);
            /* calculate the distances from node3 to other node in the Set except node1 and node2 */
            N node3 = makeNewNodeInGraph(resultTree);
            for(N a: Node){
                D ave = divide(add(getDistance(a, node1), getDistance(a, node2)), makeD(2));
                setDistance(a, node3, ave);
                setDistance(node3, a, ave);
            }
            /* add new node to Set node */
            Node.add(node3);
            /* add new edge to Graph */
            addEdgeToGraph(node1, node3, resultTree);
            addEdgeToGraph(node2, node3, resultTree);
        }

       return resultTree;
    }

    /***
     * Returns the distance between two distinct taxa.
     *
     * @param taxon1 a taxon
     * @param taxon2 another taxon
     * @param dist
     * @return the distance between teh two taxa.
     */
    protected abstract D setDistance(N taxon1, N taxon2, D dist);


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
