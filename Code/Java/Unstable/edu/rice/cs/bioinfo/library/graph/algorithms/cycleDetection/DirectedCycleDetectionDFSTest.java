package edu.rice.cs.bioinfo.library.graph.algorithms.cycleDetection;

import junit.framework.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/9/13
 * Time: 5:45 PM
 * To change this template use File | Settings | File Templates.
 */
public class DirectedCycleDetectionDFSTest
{
    @Test
    public void testTryFindCycle()
    {
        Set<String> edges = new HashSet<String>(Arrays.asList("AB", "BC", "CA"));
        List<String> cycle = new DirectedCycleDetectionDFS<String>()
        {
            @Override
            protected Object getDestination(String edge) {
                 return edge.charAt(1);
            }

            @Override
            protected Object getSource(String edge) {
                return edge.charAt(0);
            }
        }.tryFindCycle(edges);

        Assert.assertEquals(3, cycle.size());
        Assert.assertTrue(cycle.contains("AB"));
        Assert.assertTrue(cycle.contains("BC"));
        Assert.assertTrue(cycle.contains("CA"));

    }
}
