package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/15/13
 * Time: 5:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class FindAllMinimumSpanningTreesCombinationsAscTest extends FindAllMinimumSpanningTreesTest
{
    class FindAllMinimumSpanningTreesCombinationsAscImpl extends
            FindAllMinimumSpanningTreesCombinationsAsc<Set<Tuple<String, String>>,String,
                    Tuple<String,String>, Integer>

    {
        private Set<Tuple<String, String>> _edges;

        private Map<Tuple<String,String>,Integer> _edgeToWeight;

        FindAllMinimumSpanningTreesCombinationsAscImpl(Set<Tuple<String, String>> edges,
                                                       Map<Tuple<String,String>,Integer> edgeToWeight) throws GraphDisconnectedException
        {
            super(edges);
            _edges = edges;
            _edgeToWeight = edgeToWeight;
            setup();
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
            return graph;
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
        protected Integer getWeight(Tuple<String, String> edge)
        {
            return _edgeToWeight.get(edge);
        }
    }
    /*
    @Test
    public void testExecute2Nodes1Edge() throws GraphDisconnectedException
    {
        Set<Tuple<String,String>> edges = new HashSet<Tuple<String, String>>();
        Tuple<String,String> e1 = new Tuple<String, String>("A", "B");
        edges.add(e1);

        HashMap<Tuple<String,String>,Integer> edgeToWeight = new HashMap<Tuple<String, String>, Integer>();
        edgeToWeight.put(e1, 1);

        int numMsts = 0;
        for(Set<Tuple<String,String>> mst : new FindAllMinimumSpanningTreesCombinationsAscImpl(edges, edgeToWeight))
        {
            numMsts++;
        }

        Assert.assertTrue(numMsts == 1);
    }

    @Test
    public void testExecute3NodeComplete() throws GraphDisconnectedException
    {
        Set<String> nodes = new HashSet<String>();
        nodes.add("A");
        nodes.add("B");
        nodes.add("C");

        final HashMap<Tuple<String,String>,Integer> edgeToWeight = new HashMap<Tuple<String, String>, Integer>();

        Set<Tuple<String,String>> graphEdges = new CompleteGraphFactory<String,Tuple<String,String>>()
        {
            @Override
            public Tuple<String, String> makeEdge(String node1, String node2)
            {
                Tuple<String,String> edge = new Tuple<String, String>(node1, node2);
                edgeToWeight.put(edge, 1);
                return edge;
            }
        }.makeCompleteGraph(nodes);

        int numMsts = 0;
        for(Set<Tuple<String,String>> mst : new FindAllMinimumSpanningTreesCombinationsAscImpl(graphEdges, edgeToWeight))
        {
            numMsts++;
        }

        Assert.assertTrue(numMsts == 3);
    }

    @Test
    public void testExecute4NodeComplete() throws GraphDisconnectedException
    {
        Set<String> nodes = new HashSet<String>();
        nodes.add("R");
        nodes.add("G");
        nodes.add("B");
        nodes.add("Y");

        final HashMap<Tuple<String,String>,Integer> edgeToWeight = new HashMap<Tuple<String, String>, Integer>();

        Set<Tuple<String,String>> graphEdges = new CompleteGraphFactory<String,Tuple<String,String>>()
        {
            @Override
            public Tuple<String, String> makeEdge(String node1, String node2)
            {
                Tuple<String,String> edge = new Tuple<String, String>(node1, node2);
                edgeToWeight.put(edge, 1);
                return edge;
            }
        }.makeCompleteGraph(nodes);

        int numMsts = 0;
        for(Set<Tuple<String,String>> mst : new FindAllMinimumSpanningTreesCombinationsAscImpl(graphEdges, edgeToWeight))
        {
            numMsts++;
        }

        Assert.assertTrue(numMsts == 16);
    }      */

    @Override
    protected Iterable<Set<Tuple<String, String>>> makeSearcher(Set<String> nodes, Set<Tuple<String, String>> edges, Map<Tuple<String, String>, Integer> edgeToWeight) throws GraphDisconnectedException
    {
        return new FindAllMinimumSpanningTreesCombinationsAscImpl(edges, edgeToWeight);
    }


}
