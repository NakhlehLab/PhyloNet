package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import edu.rice.cs.bioinfo.library.graph.algorithms.bridge.IsBridgeDfs;
import edu.rice.cs.bioinfo.library.graph.algorithms.cycleDetection.CycleDetectorFromEdgesSets;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/10/13
 * Time: 3:48 PM
 * To change this template use File | Settings | File Templates.
 */
abstract class BacktrackSpanningTrees<G,N,E>
{
    class IsBridge extends IsBridgeDfs<Set<E>,E,N>
    {
        private final G _graph;

        IsBridge(G graph)
        {
            _graph = graph;
        }


        @Override
        protected int getNodeCount(Set<E> graph)
        {
            Set<N> nodes = new HashSet<N>();
            for(E edge : graph)
            {
                Tuple<? extends N, ? extends N> nodesOfEdge = getNodesOfEdge(edge, graph);
                nodes.add(nodesOfEdge.Item1);
                nodes.add(nodesOfEdge.Item2);

            }

            return nodes.size();
        }

        @Override
        protected Iterable<? extends E> getIncidentEdges(N node, Set<E> graph)
        {
            Set<E> incidentEdges = new HashSet<E>();
            for(E edge : graph)
            {
                Tuple<? extends N, ? extends N> nodesOfEdge = getNodesOfEdge(edge, graph);
                if(nodesOfEdge.Item1.equals(node) || nodesOfEdge.Item2.equals(node))
                    incidentEdges.add(edge);

            }

            return incidentEdges;
        }

        @Override
        protected Tuple<? extends N, ? extends N> getNodesOfEdge(E edge, Set<E> graph)
        {
            return BacktrackSpanningTrees.this.getNodesOfEdge(edge, _graph);
        }
    }

    private Set<Proc1<Set<E>>> _spanningTreeFoundListeners = new HashSet<Proc1<Set<E>>>();


    public void addSpanningTreeFoundListener(Proc1<Set<E>> listener )
    {
        _spanningTreeFoundListeners.add(listener);
    }

    public void execute(final G graph)
    {
        Set<E> potentialSpanTreeEdges = performGetEdgesOfGraph(graph);
        Set<N> uncoveredNodes = getNodesOfGraph(graph);
        int numEdgesInSpanningTree = uncoveredNodes.size() - 1;


        Set<E> partialSpan = findAllBridges(potentialSpanTreeEdges, graph);

        coverNodes(uncoveredNodes, partialSpan, graph);
        potentialSpanTreeEdges.removeAll(partialSpan);

        Map<E,Tuple<? extends N, ? extends N>> edgeToNodesOfEdge = new HashMap<E,Tuple<? extends N, ? extends N>>();
        for(E edge : potentialSpanTreeEdges)
        {
            edgeToNodesOfEdge.put(edge, getNodesOfEdge(edge, graph));
        }

        CycleDetectorFromEdgesSets<E> cycleDetector = new CycleDetectorFromEdgesSets<E>()
        {
            @Override
            protected Tuple<?, ?> getNodesOfEdge(E edge)
            {
                return BacktrackSpanningTrees.this.getNodesOfEdge(edge, graph);
            }
        };


        executeHelp(potentialSpanTreeEdges, partialSpan, uncoveredNodes, graph, numEdgesInSpanningTree, cycleDetector, edgeToNodesOfEdge);
    }

    private void coverNodes(Set<N> uncoveredNodes, Set<E> edges, G graph)
    {
        for(E edge : edges)
        {
            Tuple<? extends N, ? extends N> nodesOfEdge = getNodesOfEdge(edge, graph);
            uncoveredNodes.remove(nodesOfEdge.Item1);
            uncoveredNodes.remove(nodesOfEdge.Item2);
        }
    }

    private Set<E> findAllBridges(Set<E> graph, G originalGraph)
    {
        IsBridge isBridge = new IsBridge(originalGraph);

        HashSet<E> bridges = new HashSet<E>();

        for(E edge : graph)
        {
            if(isBridge.execute(graph, edge))
                bridges.add(edge);
        }

        return bridges;
    }

    private void executeHelp(Set<E> potentialSpanTreeEdges, Set<E> partialSpan, Set<N> uncoveredNodes, G graph,
                             int numEdgesInSpanningTree, CycleDetectorFromEdgesSets<E> cycleDetector,  Map<E,Tuple<? extends N, ? extends N>> edgeToNodesOfEdge)
    {
        if(partialSpan.size() > numEdgesInSpanningTree)
            return;

        if(uncoveredNodes.isEmpty() && partialSpan.size() == numEdgesInSpanningTree)
        {
            if(!cycleDetector.containsCycle(partialSpan))
            {
                reportSpanTreeFound(new HashSet<E>(partialSpan));
                return;
            }
        }

        if(!exploreSubtree(partialSpan))
        {
            return;
        }

        if(potentialSpanTreeEdges.isEmpty())
        {
            return;
        }

        Iterator<E> potentialSpanTreeEdgesElements = potentialSpanTreeEdges.iterator();
        E edge = potentialSpanTreeEdgesElements.next();
        potentialSpanTreeEdges.remove(edge);

        // -- Assume edge included

        partialSpan.add(edge);
        Tuple<? extends N, ? extends N> nodesOfIncludedEdge = edgeToNodesOfEdge.get(edge);
        Set<N> removedNodesFromUncoveredNodes = new HashSet<N>();

        if(uncoveredNodes.remove(nodesOfIncludedEdge.Item1))
            removedNodesFromUncoveredNodes.add(nodesOfIncludedEdge.Item1);

        if(uncoveredNodes.remove(nodesOfIncludedEdge.Item2))
            removedNodesFromUncoveredNodes.add(nodesOfIncludedEdge.Item2);

        Set<E> cycleEdges = new HashSet<E>();

        potentialSpanTreeEdgesElements = potentialSpanTreeEdges.iterator();
        while(potentialSpanTreeEdgesElements.hasNext())
        {
            E potentialEdge = potentialSpanTreeEdgesElements.next();
            partialSpan.add(potentialEdge);
            if(cycleDetector.containsCycle(partialSpan))
            {
                cycleEdges.add(potentialEdge);
                potentialSpanTreeEdgesElements.remove();
            }
            partialSpan.remove(potentialEdge);
        }

        executeHelp(potentialSpanTreeEdges, partialSpan, uncoveredNodes, graph, numEdgesInSpanningTree, cycleDetector, edgeToNodesOfEdge);

        // -- Reset

        potentialSpanTreeEdges.addAll(cycleEdges);
        uncoveredNodes.addAll(removedNodesFromUncoveredNodes);
        partialSpan.remove(edge);


        // -- Assume edge excluded

        HashSet<E> partialAndPotentialEdges = new HashSet<E>();
        partialAndPotentialEdges.addAll(potentialSpanTreeEdges);
        partialAndPotentialEdges.addAll(partialSpan);
        Set<E> newBridgeEdges = findAllBridges(partialAndPotentialEdges, graph);
        newBridgeEdges.removeAll(partialSpan);
        Set<N> newNodesCoveredByNewBridgeEdges = new HashSet<N>();
        for(E newBridgeEdge : newBridgeEdges)
        {
            Tuple<? extends N, ? extends N> nodesOfEdge = edgeToNodesOfEdge.get(newBridgeEdge);
            if(uncoveredNodes.remove(nodesOfEdge.Item1))
                newNodesCoveredByNewBridgeEdges.add(nodesOfEdge.Item1);

            if(uncoveredNodes.remove(nodesOfEdge.Item2))
                newNodesCoveredByNewBridgeEdges.add(nodesOfEdge.Item2);
        }
        partialSpan.addAll(newBridgeEdges);
        potentialSpanTreeEdges.removeAll(newBridgeEdges);

        executeHelp(potentialSpanTreeEdges, partialSpan, uncoveredNodes, graph, numEdgesInSpanningTree, cycleDetector, edgeToNodesOfEdge);

        // -- Reset

        potentialSpanTreeEdges.addAll(newBridgeEdges);
        partialSpan.removeAll(newBridgeEdges);
        uncoveredNodes.addAll(newNodesCoveredByNewBridgeEdges);


        // -- Clean up

        potentialSpanTreeEdges.add(edge);
    }

    protected boolean exploreSubtree(Set<E> partialSpanTree)
    {
        return true;
    }

    private void reportSpanTreeFound(Set<E> spanTree)
    {
        for(Proc1<Set<E>> listener : _spanningTreeFoundListeners)
        {
            listener.execute(spanTree);
        }
    }



    protected abstract Set<N> getNodesOfGraph(G graph);

    protected Set<E> performGetEdgesOfGraph(G graph)
    {
        return getEdgesOfGraph(graph);
    }

    protected abstract Set<E> getEdgesOfGraph(G graph);

    protected abstract Iterable<? extends E> getIncidentEdges(N node, G graph);

    protected abstract Tuple<? extends N, ? extends N> getNodesOfEdge(E edge, G graph);
}
