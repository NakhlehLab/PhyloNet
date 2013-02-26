package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func2;


/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/10/12
 * Time: 11:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class GetOutDegree<N,E> implements Func2<GraphReadOnly<N,E>,N,Integer>
{
    public Integer execute(GraphReadOnly<N,E> graph, N node)
    {
        if(!graph.isRooted())
        {
            throw new IllegalArgumentException("Given graph must be rooted");
        }

        int outDegree = 0;
        for(E edge : graph.getIncidentEdges(node))
        {
            if(graph.getNodesOfEdge(edge).Item1.equals(node))
            {
                outDegree++;
            }
        }

        return outDegree;
    }
}
