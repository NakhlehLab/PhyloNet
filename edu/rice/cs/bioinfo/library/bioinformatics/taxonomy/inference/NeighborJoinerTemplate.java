package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference;

import java.util.*;

public abstract class NeighborJoinerTemplate<N,E,G,D extends Comparable<D>> implements NeighborJoiner<N,G>
{
    public G performJoin(Set<N> taxa)
    {
        G resultTree = makeEdgelessGraph(taxa);
        /* create a linked hash set to store the nodes */
        Set<N> node = new HashSet<N>();
        node.addAll(taxa);
        /* create a edge-distance map to store the weights between nodes */
        Map<N,Map<N,D>> nodeToNodeDistances = new HashMap<N, Map<N, D>>();
        for(N taxon1 : taxa)
        {
            Map<N,D> taxon1NodeDistances = new HashMap<N, D>();
            nodeToNodeDistances.put(taxon1, taxon1NodeDistances);
            for(N taxon2 : taxa)
            {
                if(taxon1 != taxon2){
                    D distanceFromTaxon1ToTaxon2 = getDistance(taxon1, taxon2);
                    taxon1NodeDistances.put(taxon2, distanceFromTaxon1ToTaxon2);
                }
            }
        }
        while(node.size() > 1 ) {
            /* find the minimum distance and the corresponding nodes from the edge-node map */
            Iterator<N> it = node.iterator();
            N node1 = it.next();
            N node2 = it.next();
            D temp = nodeToNodeDistances.get(node1).get(node2);

            for(N a: node){
                for(N b: node){
                    if(a != b){
                        D distab = nodeToNodeDistances.get(a).get(b);
                        if(temp.compareTo(distab) > 0){
                            temp = distab;
                            node1 = a;
                            node2 = b;
                        }
                    }
                }
            }
            /* delete node1 and node2 from Set */
            node.remove(node1);
            node.remove(node2);
           /* calculate and reset the distances from node3 to other node in the Set except node1 and node2 */
            N node3 = makeNewNodeInGraph(resultTree);
            Map<N,D> node3NodeDistances = new HashMap<N, D>();
            nodeToNodeDistances.put(node3, node3NodeDistances);
            for(N a: node){
                D ave = divide(add(nodeToNodeDistances.get(node1).get(a), nodeToNodeDistances.get(node2).get(a)), makeD(2));
                nodeToNodeDistances.get(node3).put(a, ave);
                nodeToNodeDistances.get(a).put(node3, ave);
            }
            /* delete the distances related to node1 and node2 from map */
            nodeToNodeDistances.remove(node1);
            nodeToNodeDistances.remove(node2);
            for(N a: node){
                nodeToNodeDistances.get(a).remove(node1);
                nodeToNodeDistances.get(a).remove(node2);
            }
            /* add new node to Set node */
            node.add(node3);
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
     * @return the distance between teh two taxa.
     */
    protected abstract D getDistance(N taxon1, N taxon2);

    /**
     * Makes a new graph of type G containing the given nodes.
     *
     * Returned graph will have no edges.
     *
     * @return an an edgeless graph containing the given nodes.
     */
    protected abstract G makeEdgelessGraph(Set<N> nodes);

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
     * @return the substract of the given values
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
