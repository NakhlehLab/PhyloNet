package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.edmondsAlgo;

import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.generation.AllBranchingsGeneratorBruteForce;
import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.AllDiGraphsGenerator;
import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.CompleteDigraphFactory;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import junit.framework.Assert;
import org.junit.Test;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/7/13
 * Time: 7:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class EdmondsAlgoGibbons85Test
{
    class EdmondsOverWeights extends EdmondsAlgoGibbons85<Character,String,Integer>
    {
        int nextVertex = 0;

        private final Map<String,Integer> _edgeToWeight;

        public EdmondsOverWeights(Map<String,Integer> edgeToWeight) {
            super(0,1);
            _edgeToWeight = edgeToWeight;
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
        protected String makeEdge(Character source, Character destination) {
            return source.toString() + destination.toString();
        }

        @Override
        protected Character makeVertex() {
            return ((nextVertex++) + "").charAt(0);
        }

        @Override
        protected Integer getWeightOfEdge(String edge) {
            return _edgeToWeight.get(edge);
        }

        @Override
        protected Character getDestination(String edge) {
            return edge.charAt(1);  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        protected Character getSource(String edge) {
            return edge.charAt(0);  //To change body of implemented methods use File | Settings | File Templates.
        }

        public BranchingResult invokeFindAMaximumBranching(Set<Character> nodes)
        {
            return this.findAMaximumBranching(nodes, _edgeToWeight.keySet());
        }

        public BranchingResult invokeFindAMinimumSpanningTreeOfCompleteGraph()
        {
            return this.findAMinimumSpanningTreeOfCompleteGraph(_edgeToWeight.keySet());
        }
    }

    @Test
    public void testFindAMaximumBranching1()
    {
        Set<Character> nodes = new HashSet<Character>(Arrays.asList('A', 'B'));
        Map<String,Integer> edgeToWeight = new HashMap<String, Integer>();
        edgeToWeight.put("AB", 2);
        edgeToWeight.put("BA", 1);

        EdmondsAlgoGibbons85.BranchingResult br = new EdmondsOverWeights(edgeToWeight).invokeFindAMaximumBranching(nodes);
        Assert.assertEquals(2, br.BranchWeight);
        Assert.assertEquals(1, br.Branching.size());
        Assert.assertEquals(true, br.Branching.contains("AB"));
    }

    @Test
    public void testFindAMaximumBranching2()
    {
        Set<Character> nodes = new HashSet<Character>(Arrays.asList('A', 'B'));
        Map<String,Integer> edgeToWeight = new HashMap<String, Integer>();
        edgeToWeight.put("AB", 1);
        edgeToWeight.put("BA", 2);

        EdmondsAlgoGibbons85.BranchingResult br = new EdmondsOverWeights(edgeToWeight).invokeFindAMaximumBranching(nodes);
        Assert.assertEquals(2, br.BranchWeight);
        Assert.assertEquals(1, br.Branching.size());
        Assert.assertEquals(true, br.Branching.contains("BA"));
    }

    @Test
    public void testFindAMaximumBranching3()
    {
        Set<Character> nodes = new HashSet<Character>(Arrays.asList('A', 'B', 'C'));
        Map<String,Integer> edgeToWeight = new HashMap<String, Integer>();
        edgeToWeight.put("BA", 62);
        edgeToWeight.put("CA", 54);
        edgeToWeight.put("BC", 69);
        edgeToWeight.put("CB", 44);

        EdmondsAlgoGibbons85.BranchingResult br = new EdmondsOverWeights(edgeToWeight).invokeFindAMaximumBranching(nodes);
        Assert.assertEquals(131, br.BranchWeight);
        Assert.assertEquals(2, br.Branching.size());
        Assert.assertEquals(true, br.Branching.contains("BC"));
        Assert.assertEquals(true, br.Branching.contains("BA"));
    }


    @Test
    public void testFindAMaximumBranching4()
    {
        Set<Character> nodes = new HashSet<Character>(Arrays.asList('A', 'B', 'C', 'D', 'E', 'F', 'H'));
        Map<String,Integer> edgeToWeight = new HashMap<String, Integer>();
        edgeToWeight.put("AB", 2);
        edgeToWeight.put("BC", 3);
        edgeToWeight.put("CA", 4);
        edgeToWeight.put("HC", 1);
        edgeToWeight.put("CE", 4);
        edgeToWeight.put("DA", 1);
        edgeToWeight.put("DE", 5);
        edgeToWeight.put("EF", 4);
        edgeToWeight.put("FD", 2);


        EdmondsAlgoGibbons85.BranchingResult br = new EdmondsOverWeights(edgeToWeight).invokeFindAMaximumBranching(nodes);
        Assert.assertEquals(17, br.BranchWeight);
        Assert.assertEquals(5, br.Branching.size());
        Assert.assertEquals(true, br.Branching.contains("CE"));
        Assert.assertEquals(true, br.Branching.contains("EF"));
        Assert.assertEquals(true, br.Branching.contains("FD"));
        Assert.assertEquals(true, br.Branching.contains("CA"));
        Assert.assertEquals(true, br.Branching.contains("BC"));
    }

    @Test
    public void testFindAMinimumSpanningTreeOfCompleteGraph1()
    {
        Set<Character> nodes = new HashSet<Character>(Arrays.asList('A', 'B'));
        Map<String,Integer> edgeToWeight = new HashMap<String, Integer>();
        edgeToWeight.put("AB", 1);
        edgeToWeight.put("BA", 2);

        EdmondsAlgoGibbons85.BranchingResult br = new EdmondsOverWeights(edgeToWeight).invokeFindAMinimumSpanningTreeOfCompleteGraph();
        Assert.assertEquals(1, br.BranchWeight);
        Assert.assertEquals(1, br.Branching.size());
        Assert.assertEquals(true, br.Branching.contains("AB"));
    }


    @Test
    public void testFindAMaximumBranchingDynamic()
    {
        Random rand = new Random(13);
        int numTrialsPerformed = 0;
        int trialsBudget = 70000;

        for(int numNodes = 1; true; numNodes++)
        {
            HashSet<Integer> digraphNodes = new HashSet<Integer>();
            for(int i = 0; i<numNodes; i++)
            {
                digraphNodes.add(i);
            }

            AllDiGraphsGenerator<Integer> digraphs = new AllDiGraphsGenerator<Integer>(digraphNodes);

            for(Set<Tuple<Integer,Integer>> diGraphEdges : digraphs)
            {
                testBranchingDynamicHelp(true, digraphNodes, diGraphEdges, rand);
                numTrialsPerformed++;

                if(numTrialsPerformed == trialsBudget)
                    return;
            }
        }
    }

    @Test
    public void findAMinimumSpanningTreeOfCompleteGraph()
    {
        Random rand = new Random(13);
        int numTrialsPerformed = 0;
        int trialsBudget = 7;

        for(int numNodes = 1; true; numNodes++)
        {
            HashSet<Integer> completeDigraphNodes = new HashSet<Integer>();
            for(int i = 0; i<numNodes; i++)
            {
                completeDigraphNodes.add(i);
            }

            Set<Tuple<Integer,Integer>> completeDigraphEdges = new CompleteDigraphFactory<Integer>().makeCompleteDigraph(completeDigraphNodes);

            testBranchingDynamicHelp(true, completeDigraphNodes, completeDigraphEdges, rand);
            numTrialsPerformed++;

            if(numTrialsPerformed == trialsBudget)
                return;

        }
    }


    private void testBranchingDynamicHelp(boolean maximize, Set<Integer> digraphNodes, Set<Tuple<Integer,Integer>> digraphEdges, Random rand)
    {
        final Map<Tuple<Integer,Integer>,Integer> edgeToWeight = new HashMap<Tuple<Integer, Integer>, Integer>();
        for(Tuple<Integer,Integer> edge : digraphEdges)
        {
            edgeToWeight.put(edge, rand.nextInt(100));
        }

        AllBranchingsGeneratorBruteForce<Tuple<Integer,Integer>> allBranchings = new AllBranchingsGeneratorBruteForce<Tuple<Integer, Integer>>(digraphEdges) {
            @Override
            protected Object getSource(Tuple<Integer, Integer> edge) {
                return edge.Item1;
            }

            @Override
            protected Object getDestination(Tuple<Integer, Integer> edge) {
                return edge.Item2;
            }
        };

        int maxSeenBranchingWeight = Integer.MIN_VALUE;
        int minSeenSpanningBranchingWeight = Integer.MAX_VALUE;
        for(Set<Tuple<Integer,Integer>> branching : allBranchings)
        {
            int branchWeight = 0;
            for(Tuple<Integer,Integer> edgeOfBranching : branching)
            {
                branchWeight += edgeToWeight.get(edgeOfBranching);
            }
            if(branchWeight > maxSeenBranchingWeight)
            {
                maxSeenBranchingWeight = branchWeight;
            }
            if(branchWeight < minSeenSpanningBranchingWeight && isSpanningTree(branching, digraphNodes))
            {
                minSeenSpanningBranchingWeight = branchWeight;
            }
        }

        EdmondsAlgoGibbons85<Integer,Tuple<Integer,Integer>,Integer> edmondsSolver = new EdmondsAlgoGibbons85<Integer,Tuple<Integer,Integer>,Integer>(0,1)
        {

            @Override
            protected Integer addWeight(Integer w1, Integer w2) {
                return w1 + w2;
            }

            @Override
            protected Integer subtractWeight(Integer w1, Integer w2) {
                return w1 - w2;
            }

            @Override
            protected Tuple<Integer, Integer> makeEdge(Integer source, Integer destination) {
                return new Tuple<Integer, Integer>(source, destination);
            }

            private int _nextVertexNumber = -1;
            @Override
            protected Integer makeVertex() {
                return _nextVertexNumber--;
            }

            @Override
            protected Integer getWeightOfEdge(Tuple<Integer, Integer> edge) {
                return edgeToWeight.get(edge);
            }

            @Override
            protected Integer getDestination(Tuple<Integer, Integer> edge) {
                return edge.Item2;
            }

            @Override
            protected Integer getSource(Tuple<Integer, Integer> edge) {
                return edge.Item1;
            }
        };


        if(maximize)
        {
            EdmondsAlgoGibbons85<Integer,Tuple<Integer,Integer>,Integer>.BranchingResult edmondsResult = edmondsSolver.findAMaximumBranching(digraphNodes, digraphEdges);
            Assert.assertEquals(maxSeenBranchingWeight, edmondsResult.BranchWeight.intValue());
        }
        else
        {
            EdmondsAlgoGibbons85<Integer,Tuple<Integer,Integer>,Integer>.BranchingResult edmondsResult = edmondsSolver.findAMinimumSpanningTreeOfCompleteGraph(digraphEdges);
            if(minSeenSpanningBranchingWeight != edmondsResult.BranchWeight.intValue())
            {
                edmondsSolver.findAMinimumSpanningTreeOfCompleteGraph(digraphEdges);
            }
            if(!isSpanningTree(edmondsResult.Branching, digraphNodes))
            {
                edmondsSolver.findAMinimumSpanningTreeOfCompleteGraph(digraphEdges);
            }
            Assert.assertEquals(minSeenSpanningBranchingWeight, edmondsResult.BranchWeight.intValue());
            Assert.assertTrue(isSpanningTree(edmondsResult.Branching, digraphNodes));
        }
    }

    private <T> boolean isSpanningTree(Set<Tuple<T, T>> edges, Set<T> nodes)
    {
       Set<T> possibleRoots = new HashSet<T>(nodes);
       Map<T,Integer> nodeToNumInEdges = new HashMap<T, Integer>();

       for(T node : nodes)
       {
           nodeToNumInEdges.put(node, 0);
       }

       for(Tuple<T,T> edge : edges)
       {
           possibleRoots.remove(edge.Item2);
           nodeToNumInEdges.put(edge.Item2, nodeToNumInEdges.get(edge.Item2) + 1);

           if(nodeToNumInEdges.get(edge.Item2) > 1)
               return false;
       }

       if(possibleRoots.size() != 1)
        return false;

        return edges.size() == nodes.size() -1;


    }

}
