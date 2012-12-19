package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Tuple;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/5/12
 * Time: 5:57 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class GraphBase<N, E> implements Graph<N,E>
{
    public void addNodes(N... nodes)
    {
        for(N node : nodes)
        {
            addNode(node);
        }
    }

    public void addEdges(E... edges)
    {
        for(E edge : edges)
        {
            addEdge(edge);
        }
    }

    public boolean containsEdge(N source, N destination)
    {
        boolean isRooted = this.isRooted();
        for(E e : getEdges())
        {
            Tuple<N,N> nodesOfE = getNodesOfEdge(e);

            if(source.equals(nodesOfE.Item1) && destination.equals(nodesOfE.Item2))
            {
                return true;
            }
            else if(!isRooted && source.equals(nodesOfE.Item2) && destination.equals(nodesOfE.Item1))
            {
                return true;
            }
        }

        return false;
    }
}
