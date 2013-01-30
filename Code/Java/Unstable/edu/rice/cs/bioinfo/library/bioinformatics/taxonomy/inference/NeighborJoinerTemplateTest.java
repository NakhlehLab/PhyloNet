package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import junit.framework.Assert;
import org.junit.Test;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/30/13
 * Time: 11:07 AM
 * To change this template use File | Settings | File Templates.
 */
public class NeighborJoinerTemplateTest
{
    class Graph
    {
        public Set<String> Nodes = new HashSet<String>();

        public Set<Tuple<String,String>> Edges = new HashSet<Tuple<String,String>>();
    }

    abstract class NeighborJoinerTestTemplateImp extends NeighborJoinerTemplateDoubleDistance<String,Tuple<String,String>,Graph>
    {
        @Override
        protected Graph makeEmptyGraph() {
            return new Graph();
        }

        @Override
        protected String makeNewNodeInGraph(Graph graph) {
            String newNode = graph.Nodes.size() + "";
            graph.Nodes.add(newNode);
            return  newNode;
        }

        @Override
        protected Tuple<String, String> addEdgeToGraph(String node1, String node2, Graph graph) {

            if(!graph.Nodes.contains(node1) || !graph.Nodes.contains(node2))
            {
                throw new IllegalArgumentException();
            }

            Tuple<String,String> edge = new Tuple<String, String>(node1, node2);
            graph.Edges.add(edge);

            return edge;
        }

        @Override
        protected void removeEdgeFromGraph(Tuple<String, String> edge, Graph graph) {

            if(!graph.Edges.contains(edge))
            {
                throw new IllegalArgumentException();
            }

            graph.Edges.remove(edge);
        }
    }

    @Test
    public void testPerformJoin()
    {
        NeighborJoinerTestTemplateImp joiner = new NeighborJoinerTestTemplateImp()
        {
            @Override
            protected Double getDistance(String taxon1, String taxon2) {
                List<String> edge = Arrays.asList(taxon1,taxon2);

                if(edge.containsAll(Arrays.asList("A", "B")))
                {
                    return 7.0;
                }
                else if(edge.containsAll(Arrays.asList("A", "C")))
                {
                    return 11.0;
                }
                else if(edge.containsAll(Arrays.asList("A", "D")))
                {
                    return 14.0;
                }
                else if(edge.containsAll(Arrays.asList("B", "C")))
                {
                    return 6.0;
                }
                else if(edge.containsAll(Arrays.asList("B", "D")))
                {
                    return 9.0;
                }
                else if(edge.containsAll(Arrays.asList("C", "D")))
                {
                    return 7.0;
                }

                throw new IllegalArgumentException();
            }
        };

        Graph joinTree = joiner.performJoin(new HashSet<String>(Arrays.asList("A", "B", "C", "D")));
        assertLeafs(Arrays.asList("A", "B", "C", "D"), joinTree);
        String u = assertAndGetParent("A", "B", joinTree);
        String x = assertAndGetParent("D", "C", joinTree);
        Assert.assertEquals(u, assertAndGetParent("B", x, joinTree));
        Assert.assertEquals(u, assertAndGetParent("A", x, joinTree));
        Assert.assertEquals(u, assertAndGetParent("C", u, joinTree));
        Assert.assertEquals(u, assertAndGetParent("D", u, joinTree));
    }

    private String assertAndGetParent(String node1, String node2, Graph joinTree) {
        Set<String> commonIncident = findCommonIncidentNodes(node1, node2, joinTree);
        Assert.assertEquals(1, commonIncident.size());
        return commonIncident.iterator().next();
    }

    private void assertLeafs(List<String> taxa, Graph joinTree)
    {
        Map<String,Integer> nodeToIncidentCount = new HashMap<String,Integer>();

        for(String node : joinTree.Nodes)
        {
            nodeToIncidentCount.put(node, 0);
        }

        for(Tuple<String,String> edge : joinTree.Edges)
        {
            for(String edgeNode : Arrays.asList(edge.Item1, edge.Item2))
            {
                int currentCount = nodeToIncidentCount.get(edgeNode);
                currentCount++;
                nodeToIncidentCount.put(edgeNode, currentCount);
            }
        }

        for(String taxon : taxa)
        {
            Assert.assertEquals(1, nodeToIncidentCount.get(taxon).intValue());
        }

    }

    private Set<String> findCommonIncidentNodes(String node1, String node2, Graph joinTree)
    {
        Set<String> node1Incident = new HashSet<String>();
        Set<String> node2Incident = new HashSet<String>();

        for(Tuple<String,String> edge : joinTree.Edges)
        {
            if(edge.Item1.equals(node1))
            {
                node1Incident.add(edge.Item2);
            }
            else if(edge.Item2.equals(node1))
            {
                node1Incident.add(edge.Item1);
            }
            else if(edge.Item1.equals(node2))
            {
                node2Incident.add(edge.Item2);
            }
            else if(edge.Item2.equals(node2))
            {
                node2Incident.add(edge.Item1);
            }
        }

        node1Incident.retainAll(node2Incident);

        return node1Incident;
    }
}
