package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/14/12
 * Time: 11:00 AM
 * To change this template use File | Settings | File Templates.
 */
public class IsDestinationNode<G extends GraphReadOnly<N,E>,N,E> implements Func3<G,N,E,Boolean>
{
    public Boolean execute(G graph, N node, E edge) {

        if(graph.isRooted())
        {
            Tuple<N,N> nodesOfEdge = graph.getNodesOfEdge(edge);
            return node.equals(nodesOfEdge.Item2);
        }
        else
        {
            return false;
        }
    }
}
