package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/7/12
 * Time: 7:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class FindRoot<N,E> implements Func1<GraphReadOnly<N,E>, N>
{
    public N execute(GraphReadOnly<N,E> graph)
    {
        if(!graph.isRooted())
        {
            throw new IllegalArgumentException("Passed graph must be rooted.");
        }
        IsDestinationNode<GraphReadOnly<N,E>,N,E> isDestinationNode = new IsDestinationNode();
        for (N node : graph.getNodes()) // assert no node duplicates
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
                return node;
            }
        }

        throw new IllegalArgumentException("Given directed graph has no node with indegree 0.");
    }
}
