package edu.rice.cs.bioinfo.library.phylogenetics;
import edu.rice.cs.bioinfo.library.programming.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/7/12
 * Time: 4:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphValidator
{
    public static <N,E> void assertValidGraph(GraphReadOnly<N,E> graph)
    {
        HashSet<N> nodes = new HashSet<N>();

        N searchRoot = null;
        for (N node : graph.getNodes()) // assert no node duplicates
        {
            nodes.add(node);

            if(graph.isRooted())
            {
                boolean nodeAlwaysSourceNode = true;
                for(E edge : graph.getIncidentEdges(node))
                {
                    if(graph.isDestinationNode(node, edge))
                    {
                        nodeAlwaysSourceNode = false;
                    }
                }

                if(nodeAlwaysSourceNode)
                {
                    if(searchRoot == null)
                    {
                        searchRoot = node;
                    }
                    else
                    {
                        throw new IllegalArgumentException("Multiple nodes with indegree 0 in rooted graph.");
                    }
                }
            }
            else
            {
                searchRoot = node;
            }
        }

        if(nodes.size() > 0)
        {
            if(graph.isRooted() && searchRoot == null)
            {
                throw new IllegalArgumentException("Given directed graph has no node with indegree 0.");
            }

            HashSet<N> seenNodes =  assertNoCycle(graph, searchRoot);

            if (seenNodes.size() != nodes.size())
            {
                throw new IllegalArgumentException("Given graph is not connected.");
            }
        }
    }

    private static <N,E> HashSet<N> assertNoCycle(GraphReadOnly<N,E> graph, N searchRoot)
    {
        LinkedList<N> nodeNeighbors = new LinkedList<N>();

        for(E edge : graph.getIncidentEdges(searchRoot))
        {
            Tuple<N,N> nodesOfEdge = graph.getNodesOfEdge(edge);
            if(nodesOfEdge.Item1.equals(searchRoot))
            {
                nodeNeighbors.addLast(nodesOfEdge.Item2);
            }
            else
            {
                nodeNeighbors.addLast(nodesOfEdge.Item1);
            }
        }

        HashSet<N> seenNodes = new HashSet<N>();
        seenNodes.add(searchRoot);

        for(N nodeNeighbor : nodeNeighbors)
        {
            dfsExplore(graph, graph.isRooted(), searchRoot, nodeNeighbor, seenNodes);
        }

        return seenNodes;
    }

    private static <N,E> void dfsExplore(GraphReadOnly<N,E> tree, boolean isRooted, N parent, N child, HashSet<N> seenNodes)
    {
        if (seenNodes.contains(child))
        {
            throw new IllegalArgumentException("Given graph contains a cycle.");
        }
        seenNodes.add(child);

        LinkedList<N> searchNeighbors = new LinkedList<N>();

        for(E edge : tree.getIncidentEdges(child))
        {
            Tuple<N,N> nodesOfEdge = tree.getNodesOfEdge(edge);
            N otherNode = (N) nodesOfEdge.other(child);

            if(!otherNode.equals(parent))
            {
                if(isRooted && nodesOfEdge.Item2 == otherNode)
                {
                    searchNeighbors.addLast(otherNode);
                }
                else if(!isRooted)
                {
                    searchNeighbors.addLast(otherNode);
                }
            }
        }

        for(N nonParentChild : searchNeighbors)
        {
            dfsExplore(tree, isRooted, child, nonParentChild, seenNodes);
        }
    }


}
