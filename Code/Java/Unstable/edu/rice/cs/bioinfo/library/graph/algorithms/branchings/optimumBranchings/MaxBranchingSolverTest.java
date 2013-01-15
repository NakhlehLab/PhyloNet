package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings;

import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.generation.AllBranchingsGeneratorBruteForce;
import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.AllDiGraphsGenerator;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import junit.framework.Assert;
import org.junit.Test;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/14/13
 * Time: 5:34 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class MaxBranchingSolverTest
{
    protected abstract MaxBranchingSolver<Character,String,Integer> makeSolverCharacter(Map<String,Integer> edgeToWeight);

    protected abstract MaxBranchingSolver<Integer,Tuple<Integer,Integer>,Integer> makeSolverInteger(Map<Tuple<Integer,Integer>,Integer> edgeToWeight);

    @Test
    public void testFindAMaximumBranching1()
    {
        Set<Character> nodes = new HashSet<Character>(Arrays.asList('A', 'B'));
        Map<String,Integer> edgeToWeight = new HashMap<String, Integer>();
        edgeToWeight.put("AB", 2);
        edgeToWeight.put("BA", 1);

        BranchingResult<String,Integer> br = makeSolverCharacter(edgeToWeight).findAMaximumBranching(nodes, edgeToWeight.keySet());
        Assert.assertEquals(new Integer(2), br.BranchingWeight);
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

        BranchingResult<String,Integer> br = makeSolverCharacter(edgeToWeight).findAMaximumBranching(nodes, edgeToWeight.keySet());
        Assert.assertEquals(new Integer(2), br.BranchingWeight);
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

        BranchingResult<String,Integer> br = makeSolverCharacter(edgeToWeight).findAMaximumBranching(nodes, edgeToWeight.keySet());
        Assert.assertEquals(new Integer(131), br.BranchingWeight);
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


        BranchingResult<String,Integer> br = makeSolverCharacter(edgeToWeight).findAMaximumBranching(nodes, edgeToWeight.keySet());
        Assert.assertEquals(new Integer(17), br.BranchingWeight);
        Assert.assertEquals(5, br.Branching.size());
        Assert.assertEquals(true, br.Branching.contains("CE"));
        Assert.assertEquals(true, br.Branching.contains("EF"));
        Assert.assertEquals(true, br.Branching.contains("FD"));
        Assert.assertEquals(true, br.Branching.contains("CA"));
        Assert.assertEquals(true, br.Branching.contains("BC"));
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

        Set<Set<Tuple<Integer,Integer>>> maxSeenBranchings = new HashSet<Set<Tuple<Integer,Integer>>>();
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
                maxSeenBranchings.clear();
                maxSeenBranchings.add(branching);
            }
            else if(branchWeight == maxSeenBranchingWeight)
            {
                maxSeenBranchings.add(branching);
            }
         /*   if(branchWeight < minSeenSpanningBranchingWeight && isSpanningTree(branching, digraphNodes))
            {
                minSeenSpanningBranchingWeight = branchWeight;
            }    */
        }

        MaxBranchingSolver<Integer,Tuple<Integer,Integer>,Integer> solver = makeSolverInteger(edgeToWeight);


        if(maximize)
        {
            BranchingResult<Tuple<Integer,Integer>,Integer> result = solver.findAMaximumBranching(digraphNodes, digraphEdges);
            if(maxSeenBranchingWeight != result.BranchingWeight.intValue())
            {
                BranchingResult<Tuple<Integer,Integer>,Integer> result2 = solver.findAMaximumBranching(digraphNodes, digraphEdges);
                int i = 0;
            }
            Assert.assertEquals(new Integer(maxSeenBranchingWeight), result.BranchingWeight);
            Assert.assertTrue(maxSeenBranchings.contains(result.Branching));
        }
        else
        {     /*
            BranchingResult<Tuple<Integer,Integer>,Integer> edmondsResult = edmondsSolver.findAMinimumSpanningTreeOfCompleteGraph(digraphEdges);
            if(minSeenSpanningBranchingWeight != edmondsResult.BranchingWeight.intValue())
            {
                edmondsSolver.findAMinimumSpanningTreeOfCompleteGraph(digraphEdges);
            }
            if(!isSpanningTree(edmondsResult.Branching, digraphNodes))
            {
                edmondsSolver.findAMinimumSpanningTreeOfCompleteGraph(digraphEdges);
            }
            Assert.assertEquals(new Integer(minSeenSpanningBranchingWeight), edmondsResult.BranchingWeight);
            Assert.assertTrue(isSpanningTree(edmondsResult.Branching, digraphNodes));  */
        }
    }

}
