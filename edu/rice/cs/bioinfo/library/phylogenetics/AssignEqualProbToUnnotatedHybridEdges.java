package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/10/12
 * Time: 4:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class AssignEqualProbToUnnotatedHybridEdges<N,E>  {

    public void execute(GraphReadOnly<N,E> graph, Proc3<GraphReadOnly<N,E>, E, Double> setProb, Func2<GraphReadOnly<N,E>, E, Boolean> isEdgeProbUnset)
    {
        IsDestinationNode isDestNode = new IsDestinationNode();
        for(N node : graph.getNodes())
        {
            Iterable<E> edges = graph.getIncidentEdges(node);
            LinkedList<E> inEdges = new LinkedList<E>();

            for(E edge : edges)
            {
               if(isDestNode.execute(graph, node, edge))
               {
                  inEdges.add(edge);
               }
            }

            if(inEdges.size() == 2)
            {
                if(isEdgeProbUnset.execute(graph, inEdges.get(0)) && isEdgeProbUnset.execute(graph, inEdges.get(1)))
                {
                    setProb.execute(graph, inEdges.get(0), .5);
                    setProb.execute(graph, inEdges.get(1), .5);
                }
            }
        }
    }
}
