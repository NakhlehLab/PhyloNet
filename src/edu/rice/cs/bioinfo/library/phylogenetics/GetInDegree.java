package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func2;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/10/12
 * Time: 1:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class GetInDegree<N,E> implements Func2<GraphReadOnly<N,E>,N,Integer>
{
    public Integer execute(GraphReadOnly<N,E> graph, N node)
    {
        if(!graph.isRooted())
        {
            throw new IllegalArgumentException("Given graph must be rooted");
        }

        int inDegree = 0;
        for(E edge : graph.getIncidentEdges(node))
        {
            if(graph.getNodesOfEdge(edge).Item2.equals(node))
            {
                inDegree++;
            }
        }

        return inDegree;
    }
}
