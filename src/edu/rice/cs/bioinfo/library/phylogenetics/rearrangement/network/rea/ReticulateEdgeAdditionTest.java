package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import org.junit.Assert;
import org.junit.Test;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/9/12
 * Time: 1:28 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReticulateEdgeAdditionTest<G extends Graph<String, Tuple<String,String>>>
{
    private Func3<G,String,String, Tuple<String,String>> _makeEdge;

    private Func1<G,String> _makeNode;

    public ReticulateEdgeAdditionTest()
    {
        _makeEdge = new Func3<G, String, String, Tuple<String, String>>() {
            public Tuple<String, String> execute(G network, String source, String destination) {
                return makeEdge(network, source, destination);
            }
        };
        _makeNode = new Func1<G, String>() {

            private int _i = 0;
            public String execute(G network)
            {
               return makeNode(network, "_" + (++_i));
            }
        };
    }

    protected abstract G makeNetwork(String... nodes);

    protected abstract boolean containsEdge(G network, String source, String destination);

    protected abstract String makeNode(G network, String node);

    protected abstract Tuple<String,String> makeEdge(G network, String source, String destination);

    @Test
    public void testComputeRearrangementsWithoutValidationCase1()
    {
        G network = makeNetwork(new String[] { "R", "A", "B" },
                                new String[][] { new String[] { "R", "A" },
                                                 new String[] { "R", "B" }});

        final LinkedList<G> expectedNetworks = new LinkedList<G>();
        String[] neighborNodes = new String[] { "R", "A", "B", "_1", "_2" };
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "_1" },
                                                         new String[] { "R", "_2" },
                                                         new String[] { "_1", "A" },
                                                         new String[] { "_2", "B" } }));
        expectedNetworks.add(makeNetwork(neighborNodes, new String[][]{
                new String[]{"_1", "_2"},
                new String[]{"R", "_1"},
                new String[]{"R", "_2"},
                new String[]{"_2", "A"},
                new String[]{"_1", "B"}}));
        testComputeRearrangementsWithoutValidationHelp(network, expectedNetworks);

    }

    @Test
    public void testComputeRearrangementsWithoutValidationCase2()
    {
        G network = makeNetwork(new String[] { "R", "A", "1", "B", "C" },
                                new String[][] { new String[] { "R", "A" },
                                                 new String[] { "R", "1" },
                                                 new String[] { "1", "B" },
                                                 new String[] { "1", "C" } });

        final LinkedList<G> expectedNetworks = new LinkedList<G>();
        String[] neighborNodes = new String[] { "R", "A", "1", "B", "C", "_1", "_2" };
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "1" },
                                                         new String[] { "R", "_1" },
                                                         new String[] { "_1", "A" },
                                                         new String[] { "1", "C" },
                                                         new String[] { "1", "_2" },
                                                         new String[] { "_2", "B" } }));
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "1" },
                                                         new String[] { "R", "_2" },
                                                         new String[] { "_2", "A" },
                                                         new String[] { "1", "C" },
                                                         new String[] { "1", "_1" },
                                                         new String[] { "_1", "B" } }));
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "_1" },
                                                         new String[] { "R", "_2" },
                                                         new String[] { "_1", "A" },
                                                         new String[] { "_2", "1" },
                                                         new String[] { "1", "C" },
                                                         new String[] { "1", "B" } }));
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "_1" },
                                                         new String[] { "R", "_2" },
                                                         new String[] { "_2", "A" },
                                                         new String[] { "_1", "1" },
                                                         new String[] { "1", "C" },
                                                         new String[] { "1", "B" } }));
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "_1" },
                                                         new String[] { "R", "1" },
                                                         new String[] { "_1", "A" },
                                                         new String[] { "1", "_2" },
                                                         new String[] { "_2", "C" },
                                                         new String[] { "1", "B" } }));
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "_2" },
                                                         new String[] { "R", "1" },
                                                         new String[] { "_2", "A" },
                                                         new String[] { "1", "_1" },
                                                         new String[] { "_1", "C" },
                                                         new String[] { "1", "B" } }));
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "A" },
                                                         new String[] { "R", "_1" },
                                                         new String[] { "_1", "1" },
                                                         new String[] { "1", "C" },
                                                         new String[] { "1", "_2" },
                                                         new String[] { "_2", "B" } }));
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "A" },
                                                         new String[] { "R", "_1" },
                                                         new String[] { "_1", "1" },
                                                         new String[] { "1", "_2" },
                                                         new String[] { "_2", "C" },
                                                         new String[] { "1", "B" } }));
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "A" },
                                                         new String[] { "R", "1" },
                                                         new String[] { "1", "_1" },
                                                         new String[] { "1", "_2" },
                                                         new String[] { "_1", "C" },
                                                         new String[] { "_2", "B" } }));
        expectedNetworks.add(makeNetwork(neighborNodes,  new String[][] {
                                                         new String[] { "_1", "_2" },
                                                         new String[] { "R", "A" },
                                                         new String[] { "R", "1" },
                                                         new String[] { "1", "_1" },
                                                         new String[] { "1", "_2" },
                                                         new String[] { "_2", "C" },
                                                         new String[] { "_1", "B" } }));
        testComputeRearrangementsWithoutValidationHelp(network, expectedNetworks);


    }

    private void testComputeRearrangementsWithoutValidationHelp(G network, final Collection<G> expectedNetworks)
    {
         new ReticulateEdgeAdditionInPlace<G, String, Tuple<String,String>>(_makeNode, _makeEdge).computeRearrangements(network, false,
                new Func4<G, Tuple<String, String>, Tuple<String, String>, Tuple<String, String>,Boolean>() {
                    public Boolean execute(G network, Tuple<String, String> sourceEdge, Tuple<String, String> destinationEdge, Tuple<String, String> reticulateEdge)
                    {

                        G toBeRemoved = null;
                        for(G expectedNetwork : expectedNetworks)
                        {
                            if(areSameNetwork(network, expectedNetwork))
                            {
                                toBeRemoved = expectedNetwork;
                                break;
                            }
                        }

                        if(toBeRemoved != null)
                        {
                            expectedNetworks.remove(toBeRemoved);
                        }
                        else
                        {
                            Assert.fail("Generated unexpected network.");
                        }

                        return true;
                    }
                });

        Assert.assertTrue(expectedNetworks.size() == 0);
    }

    private G makeNetwork(String[] nodes, String[][] edges)
    {
        G network = makeNetwork();

        for(String node: nodes)
        {
            network.addNode(makeNode(network, node));
        }

        for(String[] edge : edges)
        {
            network.addEdge(_makeEdge.execute(network, edge[0], edge[1]));
        }

        return network;
    }

    private boolean areSameNetwork(G net1, G net2)
    {
        List<String> net1Nodes = IterableHelp.toList(net1.getNodes());
        List<String> net2Nodes  = IterableHelp.toList(net2.getNodes());

        if(net1Nodes.size() != net2Nodes.size())
        {
            return false;
        }

        List<Tuple<String,String>> net1Edges = IterableHelp.toList(net1.getEdges());
        List<Tuple<String,String>> net2Edges = IterableHelp.toList(net2.getEdges());

        if(net1Edges.size() != net2Edges.size())
        {
            return false;
        }

        return net1Nodes.containsAll(net2Nodes) && net1Edges.containsAll(net2Edges);


    }
}
