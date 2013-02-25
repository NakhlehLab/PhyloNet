package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/28/12
 * Time: 1:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphMask<N,E,G extends GraphReadOnly<N,E>> implements GraphReadOnly<N,E>
{
    public final G Graph;

    public final Iterable<N> ExcludedNodes;

    public final Iterable<E> ExcludedEdges;

    public GraphMask(G graph, Iterable<N> excludedNodes, Iterable<E> excludedEdges)
    {
        Graph = graph;

        HashSet<E> excludedEdgesBacking = new HashSet<E>(IterableHelp.<E,E>toList(excludedEdges));

        for(N node : excludedNodes)
        {
            for(E edge : graph.getIncidentEdges(node))
            {
                if(!excludedEdgesBacking.contains(edge))
                {
                    excludedEdgesBacking.add(edge);
                }
            }
        }


        ExcludedNodes = excludedNodes;
        ExcludedEdges = excludedEdgesBacking;
    }

    public Iterable<N> getNodes() {

        List<N> tbr = IterableHelp.toList(Graph.getNodes());
        tbr.removeAll(IterableHelp.toList(ExcludedNodes));
        return tbr;
    }

    public Iterable<E> getEdges() {
        List<E> tbr = IterableHelp.toList(Graph.getEdges());
        tbr.removeAll(IterableHelp.toList(ExcludedEdges));
        return tbr;
    }

    public Tuple<N, N> getNodesOfEdge(E edge) {

        for(E e : ExcludedEdges)
        {
            if(edge.equals(e))
            {
                throw new IllegalArgumentException("Given edge is a member of the excluded edges.");
            }
        }
        return Graph.getNodesOfEdge(edge);
    }

    public Iterable<E> getIncidentEdges(N node) {

        for(N n : ExcludedNodes)
        {
            if(node.equals(n))
            {
                throw new IllegalArgumentException("Given node is a member of the excluded nodes.");
            }
        }
        List<E> tbr = new LinkedList<E>(IterableHelp.<E,E>toList(Graph.getIncidentEdges(node)));
        tbr.removeAll(IterableHelp.toList(ExcludedEdges));
        return tbr;
    }

    public E getEdge(N source, N destination) {

        for(N n : ExcludedNodes)
        {
            if(source.equals(n))
            {
                throw new IllegalArgumentException("Given source node is a member of the excluded nodes.");
            }
        }

        for(N n : ExcludedNodes)
        {
            if(destination.equals(n))
            {
                throw new IllegalArgumentException("Given destination node is a member of the excluded nodes.");
            }
        }

        E edge = Graph.getEdge(source, destination);

        for(E e : ExcludedEdges)
        {
            if(edge.equals(e))
            {
                throw new IllegalArgumentException("Given edge is a member of the excluded edges.");
            }
        }

        return edge;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isRooted() {
       return Graph.isRooted();
    }

    public boolean containsEdge(N source, N destination) {
        return getEdge(source,destination) != null;
    }
}
