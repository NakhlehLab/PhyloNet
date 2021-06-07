package edu.rice.cs.bioinfo.library.graph.algorithms.minimumSpanningArborescence;

import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.generation.AllBranchingsGeneratorBruteForce;
import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.CompleteDigraphFactory;
import junit.framework.Assert;
import org.junit.Test;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/15/13
 * Time: 1:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class MinSpanArborescenceSolverCompleteDigraphTest
{
    class MinSpanArborescenceSolverCompleteDigraphImp extends MinSpanArborescenceSolverCompleteDigraph<WeightedEdge,Integer>
    {
        public MinSpanArborescenceSolverCompleteDigraphImp() {
            super(0, 1);
        }

        @Override
        protected Integer addWeight(Integer w1, Integer w2) {
            return w1 + w2;
        }

        @Override
        protected Integer subtractWeight(Integer w1, Integer w2) {
            return w1 - w2;
        }


        @Override
        protected Integer getWeightOfEdge(WeightedEdge edge) {
            return edge.Weight;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        protected Integer getSource(WeightedEdge edge) {
            return edge.Source;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        protected Integer getDestination(WeightedEdge edge) {
            return edge.Destination;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    class WeightedEdge
    {
        public final Integer Source;

        public final Integer Destination;

        public final Integer Weight;

        WeightedEdge(Integer source, Integer destination, Integer weight) {
            Source = source;
            Destination = destination;
            Weight = weight;
        }
    }

    @Test
    public void testTryFindMinSpanArborescence1()
    {


        List<WeightedEdge> edges = Arrays.asList(new WeightedEdge(0, 1, 1), new WeightedEdge(1, 0, 2), new WeightedEdge(1, 2, 4), new WeightedEdge(2, 1, 3), new WeightedEdge(0, 2, 5),
                new WeightedEdge(2, 0, 6));


        MinSpanArborescence<WeightedEdge,Integer> minSpanTree = new MinSpanArborescenceSolverCompleteDigraphImp().tryFindMinSpanArborescence(new HashSet<WeightedEdge>(edges));
        Assert.assertEquals(5, minSpanTree.SpanWeight.intValue());
        Assert.assertTrue((minSpanTree.Edges.contains(edges.get(0)) && minSpanTree.Edges.contains(edges.get(2)) ||
                (minSpanTree.Edges.contains(edges.get(3)) && minSpanTree.Edges.contains(edges.get(1)))));
    }

    @Test
    public void testTryFindMinSpanArborescence2()
    {

        List<WeightedEdge> edges = Arrays.asList(new WeightedEdge(0, 1, -55), new WeightedEdge(1, 0, 129), new WeightedEdge(1, 2, 13), new WeightedEdge(2, 1, -214), new WeightedEdge(0, 2, -147),
                new WeightedEdge(2, 0, 6));


        MinSpanArborescence<WeightedEdge,Integer> minSpanTree = new MinSpanArborescenceSolverCompleteDigraphImp().tryFindMinSpanArborescence(new HashSet<WeightedEdge>(edges));
        Assert.assertEquals(-361, minSpanTree.SpanWeight.intValue());
        Assert.assertTrue(minSpanTree.Edges.contains(edges.get(3)) && minSpanTree.Edges.contains(edges.get(4)));
    }

    @Test
    public void testTryFindMinSpanArborescenceDynamic()
    {
        final Random rand = new Random(45);
        int numTrialsPerformed = 0;
        int trialsBudget = 1000;
        int numTrialsPerTopology = trialsBudget / 5;

        for(int numNodes = 2; true; numNodes++)
        {
            for(int topTry = 0; topTry<numTrialsPerTopology; topTry++)
            {

                HashSet<Integer> completeDigraphNodes = new HashSet<Integer>();
                for(int i = 0; i<numNodes; i++)
                {
                    completeDigraphNodes.add(i);
                }

                Set<WeightedEdge> completeDigraphEdges = new CompleteDigraphFactory<Integer, WeightedEdge>()
                {
                    @Override
                    public WeightedEdge makeEdge(Integer source, Integer destination) {
                        return new WeightedEdge(source, destination, 50 - rand.nextInt(100));
                    }
                }.makeCompleteDigraph(completeDigraphNodes);

                AllBranchingsGeneratorBruteForce<WeightedEdge> allBranchings = new AllBranchingsGeneratorBruteForce<WeightedEdge>(completeDigraphEdges) {
                    @Override
                    protected Object getSource(WeightedEdge edge) {
                        return edge.Source;
                    }

                    @Override
                    protected Object getDestination(WeightedEdge edge) {
                        return edge.Destination;
                    }
                };

                int minSeenSpanningBranchingWeight = Integer.MAX_VALUE;
                Set<Set<WeightedEdge>> minSeenSpanTrees = new HashSet<Set<WeightedEdge>>();
                for(Set<WeightedEdge> branching : allBranchings)
                {
                    int branchWeight = 0;
                    for(WeightedEdge edgeOfBranching : branching)
                    {
                        branchWeight += edgeOfBranching.Weight;
                    }
                    if(isSpanningArborescence(branching, completeDigraphNodes))
                    {
                        if(branchWeight == minSeenSpanningBranchingWeight)
                        {
                            minSeenSpanTrees.add(branching);
                        }
                        else if(branchWeight < minSeenSpanningBranchingWeight)
                        {
                            minSeenSpanningBranchingWeight = branchWeight;
                            minSeenSpanTrees.clear();
                            minSeenSpanTrees.add(branching);
                        }
                    }
                }

                MinSpanArborescence<WeightedEdge,Integer> minSpanTree =
                        new MinSpanArborescenceSolverCompleteDigraphImp().tryFindMinSpanArborescence(new HashSet<WeightedEdge>(completeDigraphEdges));
                Assert.assertEquals(minSeenSpanningBranchingWeight, minSpanTree.SpanWeight.intValue());
                Assert.assertTrue(minSeenSpanTrees.contains(minSpanTree.Edges));



                numTrialsPerformed++;

                if(numTrialsPerformed == trialsBudget)
                    return;
            }

        }
    }

    private boolean isSpanningArborescence(Set<WeightedEdge> edges, Set<Integer> nodes)
    {
        Set<Integer> possibleRoots = new HashSet<Integer>(nodes);
        Map<Integer,Integer> nodeToNumInEdges = new HashMap<Integer, Integer>();

        for(Integer node : nodes)
        {
            nodeToNumInEdges.put(node, 0);
        }

        for(WeightedEdge edge : edges)
        {
            possibleRoots.remove(edge.Destination);
            nodeToNumInEdges.put(edge.Destination, nodeToNumInEdges.get(edge.Destination) + 1);

            if(nodeToNumInEdges.get(edge.Destination) > 1)
                return false;
        }

        if(possibleRoots.size() != 1)
            return false;

        return edges.size() == nodes.size() -1;


    }
}
