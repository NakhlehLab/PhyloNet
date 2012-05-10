package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.*;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/6/12
 * Time: 4:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class IsLeaf<N,E> implements Predicate2<GraphReadOnly<N,E>, N>
{

    public boolean execute(GraphReadOnly<N, E> tree, N node)
    {
        Iterator<E> incidentEdges = tree.getIncidentEdges(node).iterator();

        if(incidentEdges.hasNext())
        {
            E edge = incidentEdges.next();

            if(!incidentEdges.hasNext())
            {
                if(tree.isRooted() && tree.isDestinationNode(node, edge))
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
