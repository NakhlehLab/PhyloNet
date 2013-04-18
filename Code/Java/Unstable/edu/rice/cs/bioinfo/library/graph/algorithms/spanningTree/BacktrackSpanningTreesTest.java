package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import edu.rice.cs.bioinfo.library.programming.Container;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import junit.framework.Assert;
import org.junit.Test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/11/13
 * Time: 12:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class BacktrackSpanningTreesTest
{
    class BacktrackSpanningTreesTestImp extends BacktrackSpanningTrees<Map<String,Set<String>>,String,Tuple<String,String>>
    {

        @Override
        protected Set<String> getNodesOfGraph(Map<String, Set<String>> graph)
        {
            return new HashSet<String>(graph.keySet());
        }

        @Override
        protected Set<Tuple<String, String>> getEdgesOfGraph(Map<String, Set<String>> graph)
        {
            Set<Tuple<String,String>> edges = new HashSet<Tuple<String, String>>();

            for(Map.Entry<String,Set<String>> entry : graph.entrySet())
            {
                String nodeA = entry.getKey();
                for(String nodeB : entry.getValue())
                {
                    if(nodeA.compareTo(nodeB) < 0)
                        edges.add(new Tuple<String, String>(nodeA, nodeB));
                }
            }

            return  edges;

        }

        @Override
        protected Iterable<? extends Tuple<String, String>> getIncidentEdges(String node, Map<String, Set<String>> graph)
        {
            Set<Tuple<String,String>> edges = new HashSet<Tuple<String, String>>();

            for(String neighbor : graph.get(node))
            {
                edges.add(new Tuple<String, String>(node, neighbor));
            }

            return edges;
        }

        @Override
        protected Tuple<? extends String, ? extends String> getNodesOfEdge(Tuple<String, String> edge, Map<String, Set<String>> graph)
        {
            return edge;
        }
    }

    @Test
    public void testExecute1()
    {
        Map<String, Set<String>> graph = new HashMap<String, Set<String>>();

        HashSet<String> aNeighbors = new HashSet<String>();
        aNeighbors.add("B");
        graph.put("A", aNeighbors);

        HashSet<String> bNeighbors = new HashSet<String>();
        bNeighbors.add("A");
        graph.put("B", bNeighbors);

        BacktrackSpanningTreesTestImp bt = new BacktrackSpanningTreesTestImp();
        final Container<Integer> numSpans = new Container<Integer>(0);
        bt.addSpanningTreeFoundListener(new Proc1<Set<Tuple<String, String>>>()
        {
            public void execute(Set<Tuple<String, String>> spanTree)
            {
                numSpans.setContents(numSpans.getContents() + 1);
            }
        });
        bt.execute(graph);

        Assert.assertEquals(new Integer(1), numSpans.getContents());
    }

    @Test
    public void testExecute2()
    {
        Map<String, Set<String>> graph = new HashMap<String, Set<String>>();

        HashSet<String> aNeighbors = new HashSet<String>();
        aNeighbors.add("B");
        aNeighbors.add("C");
        graph.put("A", aNeighbors);

        HashSet<String> bNeighbors = new HashSet<String>();
        bNeighbors.add("A");
        graph.put("B", bNeighbors);

        HashSet<String> cNeighbors = new HashSet<String>();
        cNeighbors.add("A");
        graph.put("C", cNeighbors);

        BacktrackSpanningTreesTestImp bt = new BacktrackSpanningTreesTestImp();
        final Container<Integer> numSpans = new Container<Integer>(0);
        bt.addSpanningTreeFoundListener(new Proc1<Set<Tuple<String, String>>>()
        {
            public void execute(Set<Tuple<String, String>> spanTree)
            {
                numSpans.setContents(numSpans.getContents() + 1);
            }
        });
        bt.execute(graph);

        Assert.assertEquals(new Integer(1), numSpans.getContents());
    }

    @Test
    public void testExecute3()
    {
        Map<String, Set<String>> graph = new HashMap<String, Set<String>>();

        HashSet<String> aNeighbors = new HashSet<String>();
        aNeighbors.add("B");
        aNeighbors.add("C");
        graph.put("A", aNeighbors);

        HashSet<String> bNeighbors = new HashSet<String>();
        bNeighbors.add("A");
        bNeighbors.add("C");
        graph.put("B", bNeighbors);

        HashSet<String> cNeighbors = new HashSet<String>();
        cNeighbors.add("A");
        cNeighbors.add("B");
        graph.put("C", cNeighbors);

        BacktrackSpanningTreesTestImp bt = new BacktrackSpanningTreesTestImp();
        final Container<Integer> numSpans = new Container<Integer>(0);
        bt.addSpanningTreeFoundListener(new Proc1<Set<Tuple<String, String>>>()
        {
            public void execute(Set<Tuple<String, String>> spanTree)
            {
                numSpans.setContents(numSpans.getContents() + 1);
            }
        });
        bt.execute(graph);

        Assert.assertEquals(new Integer(3), numSpans.getContents());
    }
}
