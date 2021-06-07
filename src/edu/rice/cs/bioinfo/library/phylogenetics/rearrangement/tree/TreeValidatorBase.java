package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree;

import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphValidator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/4/12
 * Time: 3:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class TreeValidatorBase<N,E> implements TreeValidator<N,E>
{
    public void assertValidTree(GraphReadOnly<N, E> tree) {
         assertValidTreeHelp(tree);
    }

    public static <N,E> void assertValidTreeHelp(GraphReadOnly<N,E> tree)
    {
        new GraphValidator().assertValidGraph(tree);

        if(tree.isRooted())
        {
            for(N node : tree.getNodes())
            {
                int numInEdges = 0;
                for(E edge : tree.getIncidentEdges(node))
                {
                    if(tree.getNodesOfEdge(edge).Item2.equals(node))
                    {
                        numInEdges++;
                    }
                }

                if(numInEdges > 1)
                {
                    throw new IllegalArgumentException("Node: '" + node + "' has more than one parent");
                }
            }
        }


    }

}
