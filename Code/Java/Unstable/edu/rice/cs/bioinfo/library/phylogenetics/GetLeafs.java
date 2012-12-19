package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.*;
//import sun.plugin.javascript.navig4.Link;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/9/12
 * Time: 7:52 PM
 * To change this template use File | Settings | File Templates.
 */
public class GetLeafs<N> implements Func1<GraphReadOnly<N,?>, Iterable<N>>
{
    public Iterable<N> execute(GraphReadOnly<N, ?> graph) {

        IsLeaf isLeaf = new IsLeaf();

        LinkedList<N>  leafs = new LinkedList<N>();

        for(N node : graph.getNodes())
        {
            if(isLeaf.execute(graph, node))
            {
               leafs.add(node);
            }
        }

        return leafs;

    }
}
