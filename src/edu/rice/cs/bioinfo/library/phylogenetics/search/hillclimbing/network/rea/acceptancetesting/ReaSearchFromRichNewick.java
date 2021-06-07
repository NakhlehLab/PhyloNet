package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea.acceptancetesting;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphValidator;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAddition;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAdditionInPlace;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea.ReaHillClimberSteepestAscent;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import junit.framework.Assert;
import org.junit.Test;

import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 2:57 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReaSearchFromRichNewick<G extends Graph<String,PhyloEdge<String>>>
{
    protected abstract G makeNetwork(String richNewick);

    protected Func1<G, String> makeNode = new Func1<G, String>()
    {
        private int counter = 0;

        public String execute(G input1) {
            return IterableHelp.countInt(input1.getNodes()) + "";
        }
    };

    protected Func3<G, String, String, PhyloEdge<String>> makeEdge = new Func3<G, String, String, PhyloEdge<String>>()
    {
        public PhyloEdge<String> execute(G arg1, String source, String destination) {
            return new PhyloEdge<String>(source, destination);
        }
    };

    @Test
    public void testInPlace1()
    {
        G network = makeNetwork("(A,(C,B)I)R;");
        ReticulateEdgeAddition<G,String,PhyloEdge<String>> reaStrategy = new ReticulateEdgeAdditionInPlace<G, String, PhyloEdge<String>>(makeNode, makeEdge);
        ReaHillClimberSteepestAscent<G,String,PhyloEdge<String>, Double> searcher = new ReaHillClimberSteepestAscent<G, String, PhyloEdge<String>, Double>(reaStrategy, false);

        LinkedList<G> expectedGen1Neighbors = new LinkedList<G>();

        expectedGen1Neighbors.add(makeNetwork("((A,6#1)5,((C)6#1,B)I)R;"));
        expectedGen1Neighbors.add(makeNetwork("((A)6#1,((C, 6#1)5,B)I)R;"));
        expectedGen1Neighbors.add(makeNetwork("((A,6#1)5,((C,B)I)6#1)R;"));
        expectedGen1Neighbors.add(makeNetwork("((A)6#1,((C,B)I,6#1)5)R;"));
        expectedGen1Neighbors.add(makeNetwork("(A,(6#1,(C,(B)6#1)I)5)R;")); // best
        expectedGen1Neighbors.add(makeNetwork("((A,6#1)5,(C,(B)6#1)I)R;"));
        expectedGen1Neighbors.add(makeNetwork("(A,(6#1,((C)6#1,B)I)5)R;"));
        expectedGen1Neighbors.add(makeNetwork("((A)6#1,(C,(B,6#1)5)I)R;"));
        expectedGen1Neighbors.add(makeNetwork("(A,((C,6#1)5,(B)6#1)I)R;"));
        expectedGen1Neighbors.add(makeNetwork("(A,((B,6#1)5,(C)6#1)I)R;"));
        final G gen1Best = expectedGen1Neighbors.get(4);

        LinkedList<G> expectedGen2Neighbors = new LinkedList<G>();
        expectedGen2Neighbors.add(makeNetwork("((A,8#2)7,((6#1,((B)6#1,C)I)5)8#2)R;"));
        expectedGen2Neighbors.add(makeNetwork("((A,8#2)7,(6#1,(((B)6#1,C)I)8#2)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("((A,8#2)7,((C,6#1)I,((B)6#1)8#2)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("((A,8#2)7,((B)6#1,(6#1,(C)8#2)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("((A,8#2)7,(((B)8#2)6#1,(6#1,C)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("((A,8#2)7,((B)6#1,((6#1)8#2,C)I)5)R;"));

        expectedGen2Neighbors.add(makeNetwork("((A)8#2,((6#1,((B)6#1,C)I)5,8#2)7)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((B)6#1,((C,6#1)I)8#2)5,8#2)7)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,((((B)8#2)6#1,(C,6#1)I)5,8#2)7)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((B)6#1,((C)8#2,6#1)I)5,8#2)7)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((B)6#1,(C,(6#1)8#2)I)5,8#2)7)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,((((B)6#1)8#2,(C,6#1)I)5,8#2)7)R;"));
        expectedGen2Neighbors.add(makeNetwork("((A)8#2,((6#1,C)I,((B)6#1,8#2)7)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((6#1,C)I)8#2,((B)6#1,8#2)7)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((6#1)8#2,C)I,((B)6#1,8#2)7)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,((6#1,C)I,(((B)8#2)6#1,8#2)7)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,((6#1,(C)8#2)I,((B)6#1,8#2)7)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("((A)8#2,((B)6#1,((C,6#1)I,8#2)7)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((B)6#1)8#2,((C,6#1)I,8#2)7)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,((B)6#1,((C,(6#1)8#2)I,8#2)7)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((B)8#2)6#1,((C,6#1)I,8#2)7)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,((B)6#1,(((C)8#2,6#1)I,8#2)7)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("((A)8#2,((B)6#1,(C,(6#1,8#2)7)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((B)8#2)6#1,(C,(6#1,8#2)7)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,((B)6#1,((C)8#2,(6#1,8#2)7)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((B)6#1)8#2,(C,(6#1,8#2)7)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("((A)8#2,(((B,8#2)7)6#1,(C,6#1)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((B,8#2)7)6#1,((C)8#2,6#1)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("((A)8#2,((B)6#1,((C,8#2)7,6#1)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((B)8#2)6#1,((C,8#2)7,6#1)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,((B)6#1,((C,8#2)7,(6#1)8#2)I)5)R;"));
        expectedGen2Neighbors.add(makeNetwork("(A,(((B)6#1)8#2,((C,8#2)7,6#1)I)5)R;"));
        final G gen2Best = expectedGen2Neighbors.getFirst();

        final Queue<LinkedList<G>> expectedGenerations = new LinkedList<LinkedList<G>>();
        expectedGenerations.offer(expectedGen1Neighbors);
        expectedGenerations.offer(expectedGen2Neighbors);

        Func1<G, Double> getScore = new Func1<G, Double>() {

            private boolean _first = true;

            public Double execute(G input) {

                if(expectedGenerations.size() > 0)
                {
                    G toRemove = null;
                    for(G net : expectedGenerations.peek())
                    {
                        if(areSameNetwork(net, input))
                        {
                            toRemove = net;
                            break;
                        }
                    }

                    if(toRemove != null)
                    {
                        expectedGenerations.peek().remove(toRemove);
                        if(expectedGenerations.peek().size() == 0)
                        {
                            expectedGenerations.remove();
                        }
                    }
                    else if(!_first && toRemove == null)
                    {
                        Assert.fail("Unexpected network generated");
                    }
                }

                _first = false;
                if(areSameNetwork(input, gen1Best))
                {
                    return 1.0;
                }
                else if(areSameNetwork(input, gen2Best))
                {
                    return 2.0;
                }
                else
                {
                    return 0.0;
                }
            }
        };

        Comparator<Double> isBetterNetwork = new Comparator<Double>() {
            public int compare(Double o1, Double o2)
            {
                return Double.compare(o1, o2);
            }
        };

        HillClimbResult<G,Double> result = searcher.search(network, getScore, isBetterNetwork, 2);
        Assert.assertTrue(expectedGen1Neighbors.size() == 0);
        Assert.assertTrue(expectedGen2Neighbors.size() == 0);
        Assert.assertTrue(expectedGenerations.size() == 0);
        Assert.assertTrue(areSameNetwork(result.BestExaminedNetwork, gen2Best));
        Assert.assertTrue(result.BestExaminedScore == 2.0);


    }

    @Test
    public void testInPlace2()
    {
        for(int i = 0; i<1000; i++)
        {
            G network = makeNetwork("((D,(E,C)X)K,(B,A)Z)R;");
            ReticulateEdgeAddition<G,String,PhyloEdge<String>> reaStrategy = new ReticulateEdgeAdditionInPlace<G, String, PhyloEdge<String>>(makeNode, makeEdge);
            ReaHillClimberSteepestAscent<G,String,PhyloEdge<String>,Double> searcher = new ReaHillClimberSteepestAscent<G, String, PhyloEdge<String>,Double>(reaStrategy, true);



            Func1<G, Double> getScore = new Func1<G, Double>()
            {

                @Override
                public Double execute(G input) {

                     new GraphValidator().assertValidGraph(input);

                    return Double.MAX_VALUE;
                }
            };

            Comparator<Double> isBetterNetwork = new Comparator<Double>()
            {
                @Override
                public int compare(Double o1, Double o2) {
                    return 1;  //To change body of implemented methods use File | Settings | File Templates.
                }
            };

            HillClimbResult<G,Double> result = searcher.search(network, getScore, isBetterNetwork, 2);
            int j = 0;
        }


    }

    private boolean areSameNetwork(G net1, G net2)
    {
        List<String> net1Nodes = IterableHelp.toList(net1.getNodes());
        List<String> net2Nodes  = IterableHelp.toList(net2.getNodes());

        if(net1Nodes.size() != net2Nodes.size())
        {
            return false;
        }

        List net1Edges = IterableHelp.toList(net1.getEdges());
        List net2Edges = IterableHelp.toList(net2.getEdges());

        if(net1Edges.size() != net2Edges.size())
        {
            return false;
        }

        return net1Nodes.containsAll(net2Nodes) && net1Edges.containsAll(net2Edges);


    }
}
