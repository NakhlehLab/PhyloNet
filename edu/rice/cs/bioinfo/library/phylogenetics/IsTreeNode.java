package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func2;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/12
 * Time: 4:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class IsTreeNode<N,E> implements Func2<N, GraphReadOnly<N,E>,Boolean>
{
    private Func2<GraphReadOnly<N,E>,N,Integer> _getInDegreeStrategy = new GetInDegree<N, E>();

    private Func2<GraphReadOnly<N,E>,N,Integer> _getOutDegreeStrategy = new GetOutDegree<N, E>();

    public Boolean execute(N node, GraphReadOnly<N,E> network)
    {
        return _getInDegreeStrategy.execute(network, node) == 1 && _getOutDegreeStrategy.execute(network, node) >= 1;
    }
}
