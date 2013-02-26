package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Predicate2;
import edu.rice.cs.bioinfo.library.programming.counters.Counter;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/12
 * Time: 4:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class GetNetworkNodeCount<N,E,I> implements Func2<GraphReadOnly<N,E>, Counter<I>, I>
{
    private Predicate2<N, GraphReadOnly<N,E>> _isNetworkNodeStrategy = new IsNetworkNode<N, E>();

    public I execute(GraphReadOnly<N,E> network, Counter<I> counter)
    {
        counter.zero();

        for(N node : network.getNodes())
        {
            if(_isNetworkNodeStrategy.execute(node, network))
            {
                counter.increment();
            }
        }

        return counter.getCount();
    }
}
