package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Predicate2;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/12
 * Time: 4:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class IsNetworkNode<N,E> implements Predicate2<N, GraphReadOnly<N,E>>
{
    private Func2<GraphReadOnly<N,E>,N,Integer> _getInDegreeStrategy = new GetInDegree<N, E>();

    public boolean execute(N node, GraphReadOnly<N,E> network)
    {
        return _getInDegreeStrategy.execute(network, node) > 1;
    }
}
