package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference;

import edu.rice.cs.bioinfo.library.language.dot_2013_1.printing.DigraphDotPrinter;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
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

        public Set<Tuple3<String,String,Double>> Edges = new HashSet<Tuple3<String,String,Double>>();

        public String toDot()
        {
            return new DigraphDotPrinter<String, Tuple3<String,String,Double>>()
            {
                protected String getNodeLabel(String node)
                {
                    return node;
                }

                @Override
                protected String getSource(Tuple3<String, String,Double> edge) {
                    return edge.Item1;  //To change body of implemented methods use File | Settings | File Templates.
                }

                @Override
                protected String getDestination(Tuple3<String, String,Double> edge) {
                    return edge.Item2;  //To change body of implemented methods use File | Settings | File Templates.
                }


            }.toDot(Nodes, Edges);
        }

    }

    abstract class NeighborJoinerTestTemplateImp extends NeighborJoinerTemplateDoubleDistance<String,Tuple3<String,String,Double>,Graph>
    {
        @Override
        protected Graph makeEdgelessGraph(Set<String> nodes) {
            Graph g = new Graph();
            g.Nodes.addAll(nodes);
            return g;
        }

        @Override
        protected String makeNewNodeInGraph(Graph graph) {
            String newNode = graph.Nodes.size() + "";
            graph.Nodes.add(newNode);
            return  newNode;
        }

        @Override
        protected Tuple3<String, String, Double> addEdgeToGraph(String node1, String node2, Double dist, Graph graph) {

            if(!graph.Nodes.contains(node1) || !graph.Nodes.contains(node2))
            {
                throw new IllegalArgumentException();
            }

            Tuple3<String,String,Double> edge = new Tuple3<String, String, Double>(node1, node2, dist);
            graph.Edges.add(edge);

            return edge;
        }

        @Override
        protected void removeEdgeFromGraph(Tuple3<String, String, Double> edge, Graph graph) {

            if(!graph.Edges.contains(edge))
            {
                throw new IllegalArgumentException();
            }

            graph.Edges.remove(edge);
        }
    }
/*
    @Test
    public void testPerformJoin1()
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
        String dot = joinTree.toDot();
        String u = assertAndGetParent("A", "B", joinTree);
        String x = assertAndGetParent("D", "C", joinTree);
        Assert.assertEquals(u, assertAndGetParent("B", x, joinTree));
        Assert.assertEquals(u, assertAndGetParent("A", x, joinTree));
        Assert.assertEquals(x, assertAndGetParent("C", u, joinTree));
        Assert.assertEquals(x, assertAndGetParent("D", u, joinTree));
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

        for(Tuple3<String,String,Double> edge : joinTree.Edges)
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

        for(Tuple3<String,String,Double> edge : joinTree.Edges)
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
 */
    @Test
    public void testPerformJoin2()
    {
        NeighborJoinerTestTemplateImp joiner = new NeighborJoinerTestTemplateImp()
        {
            @Override
            protected Double getDistance(String taxon1, String taxon2) {
                /* the distance matrix is produced by Dingqiao's geneCalculator class */
                double[][] distance = {
                        {0, 6 , 2 ,  1 , 15 , 1, 16 , 16 , 12 , 17 , 16 , 19 , 18 , 18 , 16 , 16 , 18 , 16 , 20},
                {6 , 0 , 8 , 7 , 9 , 7 , 10 , 10 , 6 , 11 , 10 , 13 , 12 , 12 , 10 , 10 , 12 , 10 , 14},
                {2 , 8 , 0 , 1 , 17 , 1 , 18 , 18 , 14 , 19 , 18,  21 , 20 , 20 , 18 , 18 , 20 , 18 , 22},
                {1 , 7 , 1 , 0 , 16 , 0 , 17 , 17 , 13,  18,  17 , 20 , 19 , 19 , 17 , 17 , 19 , 17,  21 },
                {15 , 9 , 17 , 16 , 0 , 16 , 5 ,  5 , 15 ,  6 ,  5 ,  8  , 7 ,  7 ,  5 ,  5  , 7 ,  5  , 9},
                {1 , 7 , 1 , 0 , 16 , 0 , 17 , 17 , 13 , 18,  17 , 20 , 19 , 19 , 17 , 17 , 19 , 17 , 21},
                {16 ,10 , 18 , 17 , 5 , 17 ,  0 ,  0 , 16 ,  7 ,  0  , 3 ,  2  , 2  , 0  , 0 ,  2 ,  0 ,  4 },
                {16 , 10 , 18 , 17 , 5 , 17 ,  0 ,  0 , 16 ,  7 ,  0 ,  3 ,  2  , 2 ,  0 ,  0 ,  2 ,  0 ,  4 },
                {12 , 6 , 14 , 13 , 15 , 13 , 16 , 16 ,  0 , 17 , 16 , 19 , 18 , 18 , 16 , 16 , 18 , 16 , 20 },
                {17 , 11 , 19 , 18 , 6 , 18 ,  7 ,  7 , 17 ,  0 ,  7 , 10 ,  9 ,  9 ,  7 ,  7 ,  9 ,  7 ,  9 },
                {16 , 10 , 18 , 17 , 5 , 17 ,  0 ,  0 , 16 ,  7 ,  0 ,  3 ,  2 ,  2 ,  0 ,  0 ,  2 ,  0 ,  4 },
                {19 , 13 , 21 , 20 , 8 , 20  , 3 ,  3 , 19 , 10 ,  3 ,  0 ,  5 ,  5 ,  3 ,  3 ,  5 ,  3 ,  7 },
                {18 , 12 , 20 , 19 , 7 , 19 ,  2 ,  2 , 18 ,  9 ,  2 ,  5 ,  0 ,  0 ,  2 ,  2 ,  4 ,  2 ,  2 },
                {18 , 12 , 20 , 19 , 7 , 19 ,  2 ,  2 , 18 ,  9 ,  2 ,  5 ,  0 ,  0 ,  2 ,  2 ,  4 ,  2 ,  2 },
                {16 , 10 , 18 , 17 , 5 , 17 ,  0 ,  0 , 16 ,  7 ,  0 ,  3 ,  2 ,  2 ,  0 ,  0 ,  2 ,  0 ,  4 },
                {16 , 10 , 18 , 17 , 5 , 17 ,  0 ,  0 , 16 ,  7 ,  0 ,  3 ,  2 ,  2 ,  0 ,  0 ,  2 ,  0 ,  4 },
                {18 , 12 , 20 , 19 , 7 , 19 ,  2 ,  2 , 18 ,  9 ,  2 ,  5 ,  4 ,  4 ,  2 ,  2 ,  0 ,  2 ,  6 },
                {16 , 10 , 18 , 17 , 5 , 17 ,  0 ,  0 , 16 ,  7 ,  0 ,  3 ,  2 ,  2 ,  0 ,  0 ,  2 ,  0 ,  4 },
                {20 , 14 , 22 , 21 , 9 , 21 ,  4 ,  4 , 20 ,  9 ,  4 ,  7 ,  2 ,  2 ,  4 ,  4  , 6  , 4 ,  0 },
                };
                /* the index of patient should be the 0-18 (ordered) */
                List<String> edge = Arrays.asList(taxon1,taxon2);
                List<String> taxa = new LinkedList<String>();
                taxa.add("P1B");
                taxa.add("P1T");
                taxa.add("P2");
                taxa.add("P3");
                taxa.add("P4");
                taxa.add("P5");
                taxa.add("P6");
                taxa.add("P7");
                taxa.add("P8");
                taxa.add("P9");
                taxa.add("P10");
                taxa.add("P11");
                taxa.add("P12");
                taxa.add("P13");
                taxa.add("P14");
                taxa.add("P15");
                taxa.add("P16");
                taxa.add("P17");
                taxa.add("P18");
                for(String n1 : taxa){
                    for(String n2 : taxa){
                        if(edge.containsAll(Arrays.asList(n1, n2)) && n1 != n2 && taxon1 != taxon2){
                            Double dist = new Double(distance[taxa.indexOf(taxon1)][taxa.indexOf(taxon2)]);
                            return dist;
                        }
                    }
                }
                throw new IllegalArgumentException();
            }
        };

        Graph joinTree = joiner.performJoin(new HashSet<String>(Arrays.asList("P1B", "P1T", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18")));
        assertLeafs(Arrays.asList("P1B", "P1T", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18"), joinTree);
        String dot = joinTree.toDot();
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

        for(Tuple3<String,String,Double> edge : joinTree.Edges)
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

        for(Tuple3<String,String,Double> edge : joinTree.Edges)
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
