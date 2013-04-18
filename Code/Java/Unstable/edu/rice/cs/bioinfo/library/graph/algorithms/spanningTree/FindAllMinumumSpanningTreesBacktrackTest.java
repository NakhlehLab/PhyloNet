package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/17/13
 * Time: 1:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class FindAllMinumumSpanningTreesBacktrackTest extends FindAllMinimumSpanningTreesTest
{
    class FindAllMinumumSpanningTreesBacktrackImpl
            extends FindAllMinumumSpanningTreesBacktrack<Set<Tuple<String, String>>,String,
            Tuple<String,String>, Integer>
    {

        private Map<Tuple<String,String>,Integer> _edgeToWeight;

        FindAllMinumumSpanningTreesBacktrackImpl(Map<Tuple<String,String>,Integer> edgeToWeight) throws GraphDisconnectedException
        {
            _edgeToWeight = edgeToWeight;
        }

        @Override
        protected Set<String> getNodesOfGraph(Set<Tuple<String, String>> graph)
        {

            HashSet<String> nodes = new HashSet<String>();
            for(Tuple<String,String> edge : graph)
            {
                nodes.add(edge.Item1);
                nodes.add(edge.Item2);
            }

            return nodes;
        }

        @Override
        protected Set<Tuple<String, String>> getEdgesOfGraph(Set<Tuple<String, String>> graph)
        {
            return graph;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        protected Iterable<? extends Tuple<String, String>> getIncidentEdges(final String node, Set<Tuple<String, String>> graph)
        {
            return IterableHelp.filter(graph, new Predicate1<Tuple<String, String>>()
            {
                public boolean execute(Tuple<String, String> input)
                {
                    return input.Item1.equals(node) || input.Item2.equals(node);
                }
            });
        }

        @Override
        protected Tuple<? extends String, ? extends String> getNodesOfEdge(Tuple<String, String> edge, Set<Tuple<String, String>> graph)
        {
            return edge;
        }

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
        protected Integer getWeight(Tuple<String, String> edge1)
        {
            return _edgeToWeight.get(edge1);
        }
    }

    @Override
    protected Iterable<Set<Tuple<String, String>>> makeSearcher(Set<String> nodes, Set<Tuple<String, String>> edges, Map<Tuple<String, String>, Integer> edgeToWeight) throws GraphDisconnectedException
    {
        final Set<Set<Tuple<String, String>>> msts = new HashSet<Set<Tuple<String, String>>>();

        FindAllMinumumSpanningTreesBacktrackImpl searcher = new FindAllMinumumSpanningTreesBacktrackImpl(edgeToWeight);
        searcher.addMinSpanTreeFoundListener(new Proc1<Set<Tuple<String, String>>>()
        {
            public void execute(Set<Tuple<String, String>> input)
            {
                if(!msts.add(input))
                {
                    throw new RuntimeException("Duplicate min span tree found.");
                }
            }
        });

        searcher.execute(edges);

        return msts;

    }
}
