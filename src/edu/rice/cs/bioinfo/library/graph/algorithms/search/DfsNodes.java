package edu.rice.cs.bioinfo.library.graph.algorithms.search;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Stack;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/16/13
 * Time: 12:49 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class DfsNodes<G,N> implements Iterable<N>
{
    private final G _graph;

    private final N _searchRoot;

    DfsNodes(G graph, N searchRoot)
    {
        _graph = graph;
        _searchRoot = searchRoot;
    }

    public Iterator<N> iterator()
    {
        final HashSet<N> _seenNodes = new HashSet<N>();
        _seenNodes.add(_searchRoot);

        final Stack<N> _toExplore = new Stack();
        _toExplore.add(_searchRoot);

        return new Iterator<N>()
        {
            N _next = _searchRoot;

            boolean _hasNext = true;

            public boolean hasNext()
            {
                return _hasNext;
            }

            public N next()
            {
                N toBeReturned = _next;

                if(_toExplore.isEmpty())
                {
                    _hasNext = false;
                }
                else
                {
                    N node = _toExplore.pop();
                    _next = node;
                    for(N child : getAdjacentNodes(node))
                    {
                        boolean alreadySeen = !_seenNodes.add(child);
                        if(!alreadySeen)
                            _toExplore.push(child);
                    }
                }

                return toBeReturned;
            }

            public void remove()
            {
                throw new UnsupportedOperationException();
            }
        };
    }

    protected abstract Iterable<N> getAdjacentNodes(N node);

}
