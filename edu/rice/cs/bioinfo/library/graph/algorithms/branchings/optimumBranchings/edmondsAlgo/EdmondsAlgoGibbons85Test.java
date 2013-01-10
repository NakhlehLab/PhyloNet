package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.edmondsAlgo;

import junit.framework.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/7/13
 * Time: 7:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class EdmondsAlgoGibbons85Test
{
    abstract class EdmondsBase extends EdmondsAlgoGibbons85<Character,String,Integer>
    {
        int nextVertex = 0;

        public EdmondsBase() {
            super(0);
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
        protected Character getDestination(String edge) {
            return edge.charAt(1);  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        protected Character getSource(String edge) {
            return edge.charAt(0);  //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    class SolverForGibbonsExample extends EdmondsBase
    {
        public final Set<Character> Nodes = new HashSet<Character>(Arrays.asList('A', 'B', 'C', 'D', 'E', 'F', 'H'));

        public final Set<String> Edges = new HashSet<String>(Arrays.asList("AB", "BC", "CA", "HC", "CE", "DA", "DE", "EF", "FD"));

        @Override
        protected Integer getWeightOfEdge(String edge) {
            if(edge.equals("AB"))
                return 2;
            if(edge.equals("BC"))
                return 3;
            if(edge.equals("CA"))
                return 4;
            if(edge.equals("HC"))
                return 1;
            if(edge.equals("CE"))
                return 4;
            if(edge.equals("DA"))
                return 1;
            if(edge.equals("DE"))
                return 5;
            if(edge.equals("EF"))
                return 4;
            if(edge.equals("FD"))
                return 2;

            throw new IllegalArgumentException();
        };
    }



    @Test
    public void testTextExampleMax()
    {
        SolverForGibbonsExample gibbonsExample = new SolverForGibbonsExample();

        EdmondsAlgoGibbons85.BranchingResult br = gibbonsExample.findAMaximumBranching(gibbonsExample.Nodes, gibbonsExample.Edges);
        Assert.assertEquals(17, br.BranchWeight);
        Assert.assertEquals(5, br.Branching.size());
        Assert.assertEquals(true, br.Branching.contains("CE"));
        Assert.assertEquals(true, br.Branching.contains("EF"));
        Assert.assertEquals(true, br.Branching.contains("FD"));
        Assert.assertEquals(true, br.Branching.contains("CA"));
        Assert.assertEquals(true, br.Branching.contains("BC"));
    }

    @Test
    public void testTextExampleMin()
    {
        SolverForGibbonsExample gibbonsExample = new SolverForGibbonsExample();

        EdmondsAlgoGibbons85.BranchingResult br = gibbonsExample.findAMinimumBranching(gibbonsExample.Nodes, gibbonsExample.Edges);
        Assert.assertEquals(14, br.BranchWeight);
        Assert.assertEquals(6, br.Branching.size());
        Assert.assertEquals(true, br.Branching.contains("HC"));
        Assert.assertEquals(true, br.Branching.contains("CE"));
        Assert.assertEquals(true, br.Branching.contains("EF"));
        Assert.assertEquals(true, br.Branching.contains("FD"));
        Assert.assertEquals(true, br.Branching.contains("DA"));
        Assert.assertEquals(true, br.Branching.contains("AB"));
    }

}
