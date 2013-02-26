package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Predicate2;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/6/12
 * Time: 4:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class IsLeaf implements Predicate2<GraphReadOnly, Object>
{

    public boolean execute(GraphReadOnly tree, Object node)
    {
        Iterator incidentEdges = tree.getIncidentEdges(node).iterator();

        if(incidentEdges.hasNext())
        {
            Object edge = incidentEdges.next();

            if(!incidentEdges.hasNext())
            {
                if(tree.isRooted() && new IsDestinationNode().execute(tree, node, edge))
                {
                    return true;
                }
                else if(!tree.isRooted())
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }

    }
}
