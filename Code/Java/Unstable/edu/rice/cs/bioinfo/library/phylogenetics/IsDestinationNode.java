package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Tuple;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/14/12
 * Time: 11:00 AM
 * To change this template use File | Settings | File Templates.
 */
public class IsDestinationNode implements Func3<GraphReadOnly,Object,Object,Boolean>
{
    public Boolean execute(GraphReadOnly graph, Object node, Object edge) {

        if(graph.isRooted())
        {
            Tuple nodesOfEdge = graph.getNodesOfEdge(edge);
            return node.equals(nodesOfEdge.Item2);
        }
        else
        {
            return false;
        }
    }
}
