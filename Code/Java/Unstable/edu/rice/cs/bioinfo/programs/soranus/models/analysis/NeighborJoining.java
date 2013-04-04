package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference.NeighborJoinerTemplateIntegerDistance;
import edu.rice.cs.bioinfo.programs.soranus.models.data.Sequencing;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/5/13
 * Time: 7:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class NeighborJoining extends NeighborJoinerTemplateIntegerDistance<Sequencing, NeighborJoining.Edge, NeighborJoining.Graph>
{
    public class Edge
    {
        public final Sequencing Node1;

        public final Sequencing Node2;

        public final Integer Distance;

        public Edge(Sequencing node1, Sequencing node2, Integer distance) {
            Node1 = node1;
            Node2 = node2;
            Distance = distance;
        }
    }

    public class Graph
    {
        public final Set<Sequencing> Nodes = new HashSet<Sequencing>();

        public final Set<Edge> Edges = new HashSet<Edge>();
    }

    private Sequencing _lastNodeCreated = null;

    public Sequencing getLastNodeCreated()
    {
        return _lastNodeCreated;
    }

    @Override
    protected Integer getDistance(Sequencing taxon1, Sequencing taxon2) {
        return taxon1.getGeneticDistance(taxon2);
    }

    @Override
    protected Graph makeEdgelessGraph(Set<Sequencing> nodes) {
        Graph g = new Graph();
        g.Nodes.addAll(nodes);
        return g;
    }

    @Override
    protected Sequencing makeNewNodeInGraph(Graph graph) {
        Sequencing s = new Sequencing("");
        _lastNodeCreated = s;
        graph.Nodes.add(s);
        return s;
    }

    @Override
    protected  NeighborJoining.Edge addEdgeToGraph(Sequencing node1, Sequencing node2, Integer dist, Graph graph) {
        Edge e = new Edge(node1,node2,dist);
        graph.Edges.add(e);
        return e;

    }

    @Override
    protected void removeEdgeFromGraph(NeighborJoining.Edge edge, Graph graph) {
        graph.Edges.remove(edge);
    }
}
