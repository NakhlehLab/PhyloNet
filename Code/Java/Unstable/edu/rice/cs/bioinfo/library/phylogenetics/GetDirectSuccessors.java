package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/10/12
 * Time: 11:22 AM
 * To change this template use File | Settings | File Templates.
 */
public class GetDirectSuccessors<N,E> implements Func2<GraphReadOnly<N,E>,N,Iterable<N>>
{
    public Iterable<N> execute(GraphReadOnly<N,E> graph, N node)
    {
        if(!graph.isRooted())
        {
            throw new IllegalArgumentException("Given graph must be rooted");
        }

        LinkedList<N> directSuccessors = new LinkedList<N>();
        for(E edge : graph.getIncidentEdges(node))
        {
            Tuple<N,N> nodesOfEdge =  graph.getNodesOfEdge(edge);
            if(nodesOfEdge.Item1.equals(node))
            {
                directSuccessors.add(nodesOfEdge.Item2);
            }
        }

        return directSuccessors;
    }
}
