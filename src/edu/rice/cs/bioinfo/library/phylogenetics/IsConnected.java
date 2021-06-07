package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/28/12
 * Time: 4:34 PM
 * To change this template use File | Settings | File Templates.
 */
// weakly connected for dag
public class IsConnected  implements Predicate1<GraphReadOnly>
{
    public boolean execute(GraphReadOnly graph)
    {
        Iterator nodes = graph.getNodes().iterator();

        if(!nodes.hasNext())
        {
            return true;
        }

        Object searchRoot = nodes.next();
        Queue toExpand = new LinkedList();
        toExpand.offer(searchRoot);
        HashSet seenNodes = new HashSet();
        seenNodes.add(searchRoot);

        while(toExpand.size() > 0)
        {
           Object node = toExpand.remove();

           for(Object incidentEdge : graph.getIncidentEdges(node))
           {
               Tuple nodesOfEdge = graph.getNodesOfEdge(incidentEdge);

               for(Object neighborNode : new Object[] { nodesOfEdge.Item1, nodesOfEdge.Item2 })
               {
                    if(!seenNodes.contains(neighborNode))
                    {
                        seenNodes.add(neighborNode);
                        toExpand.add(neighborNode);
                    }
               }
           }
        }

        return seenNodes.size() == IterableHelp.countInt(graph.getNodes());
    }
}
