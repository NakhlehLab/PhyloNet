package edu.rice.cs.bioinfo.library.phylogenetics;
import edu.rice.cs.bioinfo.library.programming.*;
import org.junit.internal.matchers.Each;

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
        IsDestinationNode isDestinationNode = new IsDestinationNode();
        for (N node : graph.getNodes()) // assert no node duplicates
        {
            nodes.add(node);

            if(graph.isRooted())
            {
                boolean nodeAlwaysSourceNode = true;
                for(E edge : graph.getIncidentEdges(node))
                {
                    if(isDestinationNode.execute(graph, node, edge))
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
        HashSet<N> unfinishedNodes = new HashSet<N>();
        HashSet<N> finishedNodes = new HashSet<N>();

        boolean isGraphRooted = graph.isRooted();

        for(E incidentToRoot : graph.getIncidentEdges(searchRoot))
        {
            Tuple<N,N> nodesOfEdge = graph.getNodesOfEdge(incidentToRoot);
            N nonRoot = (N) nodesOfEdge.other(searchRoot);

            if(searchRoot.equals(nonRoot))
            {
                throw new IllegalArgumentException("Root contains a self loop.");
            }
            dfsVisit(graph, isGraphRooted, nonRoot, searchRoot, unfinishedNodes, finishedNodes);
            finishedNodes.add(searchRoot);
        }


        return finishedNodes;
    }

    private static <N,E> void dfsVisit(GraphReadOnly<N,E> graph, boolean isGraphRooted, N node, N searchParent,  HashSet<N> unfinishedNodes, HashSet<N> finishedNodes)
    {
        unfinishedNodes.add(node);

        for(E incidentEdge : graph.getIncidentEdges(node))
        {
            Tuple<N,N> nodesOfEdge = graph.getNodesOfEdge(incidentEdge);

            if(isGraphRooted && nodesOfEdge.Item2.equals(node)) // incoming edge
            {
                continue;
            }

            N otherNode =  (N) nodesOfEdge.other(node);

            if(otherNode.equals(searchParent)) // edge we just came from
            {
                continue;
            }

            if(finishedNodes.contains(otherNode)) // already fully explored this node.
            {
                continue;
            }

            if(unfinishedNodes.contains(otherNode))
            {
                throw new IllegalArgumentException("Given graph contains a cycle.");
            }

            dfsVisit(graph, isGraphRooted, otherNode, node, unfinishedNodes, finishedNodes);
        }
        finishedNodes.add(node);

    }


}
