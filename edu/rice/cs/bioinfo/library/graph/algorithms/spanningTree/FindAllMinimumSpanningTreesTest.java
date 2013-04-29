package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.CompleteGraphFactory;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import junit.framework.Assert;
import org.junit.Test;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/17/13
 * Time: 1:03 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class FindAllMinimumSpanningTreesTest
{

    protected abstract Iterable<Set<Tuple<String,String>>>
        makeSearcher(Set<String> nodes, Set<Tuple<String,String>> edges, Map<Tuple<String,String>,Integer> edgeToWeight) throws GraphDisconnectedException;

    @Test
    public void testExecuteCompleteDynamic() throws GraphDisconnectedException
    {
        for(int numNodes = 4; numNodes<=10; numNodes++)
        {
            Date start = new Date();

            Set<String> nodes = new HashSet<String>();
            for(int node = 1; node<=numNodes; node++)
            {
                nodes.add(node + "");
            }

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
            for(Set<Tuple<String,String>> mst : makeSearcher(nodes, graphEdges, edgeToWeight))
            {
                numMsts++;
            }

            Date stop = new Date();
            Date delta = new Date(stop.getTime() - start.getTime());
            long milis = delta.getTime();

            System.out.println(numNodes + " " + milis);

            double expectedNumberOfMsts = Math.pow(numNodes, numNodes-2);
            Assert.assertEquals(expectedNumberOfMsts, (double) numMsts);


        }
    }
}
