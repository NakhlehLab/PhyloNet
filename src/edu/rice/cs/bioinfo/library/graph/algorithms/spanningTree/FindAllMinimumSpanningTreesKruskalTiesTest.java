package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 5/1/13
 * Time: 3:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class FindAllMinimumSpanningTreesKruskalTiesTest extends FindAllMinimumSpanningTreesTest
{
    @Override
    protected Iterable<Set<Tuple<String, String>>> makeSearcher(final Set<String> nodes,
                                                                final Set<Tuple<String, String>> edges,
                                                                final Map<Tuple<String, String>, Integer> edgeToWeight) throws GraphDisconnectedException
    {
        return new FindAllMinimumSpanningTreesKruskalTies<Object,String,Tuple<String, String>,Integer>(nodes)
        {

            @Override
            protected Integer add(Integer term1, Integer term2)
            {
                return term1 + term2;
            }

            @Override
            protected Integer makeZero()
            {
                return 0;
            }

            @Override
            protected Set<String> getNodesOfGraph(Object graph)
            {
                return nodes;
            }

            @Override
            protected Set<Tuple<String, String>> getEdgesOfGraph(Object graph)
            {
                return edges;
            }

            @Override
            protected Iterable<? extends Tuple<String, String>> getIncidentEdges(final String node, Object graph)
            {
                return IterableHelp.filter(edges, new Predicate1<Tuple<String, String>>()
                                {
                                    public boolean execute(Tuple<String, String> input)
                                    {
                                        return input.Item1.equals(node) || input.Item2.equals(node);
                                    }
                                });
            }

            @Override
            protected Tuple<? extends String, ? extends String> getNodesOfEdge(Tuple<String, String> edge, Object graph)
            {
                return edge;
            }

            @Override
            protected Integer getWeightOfEdge(Tuple<String, String> edge)
            {
                return edgeToWeight.get(edge);
            }
        };
    }
}
