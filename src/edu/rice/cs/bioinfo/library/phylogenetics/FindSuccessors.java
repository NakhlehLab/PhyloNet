package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func2;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/18/12
 * Time: 7:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class FindSuccessors<N,E> implements Func2<GraphReadOnly<N,E>, N, Iterable<N>>
{
    @Override
    public Iterable<N> execute(GraphReadOnly<N, E> network, N node)
    {
        if(!network.isRooted())
        {
            throw new IllegalArgumentException("Given graph must be rooted");
        }

        LinkedList<N> sucessors = new LinkedList<N>();

        LinkedList<N> toSearch = new LinkedList<N>();
        toSearch.add(node);

        GetDirectSuccessors<N,E> getDirectSuccessors = new GetDirectSuccessors<N, E>();
        while(!toSearch.isEmpty())
        {
            N searchNode = toSearch.pop();

            for(N ds : getDirectSuccessors.execute(network, searchNode))
            {
                sucessors.add(ds);
                toSearch.add(ds);
            }
        }

        return sucessors;
    }
}

