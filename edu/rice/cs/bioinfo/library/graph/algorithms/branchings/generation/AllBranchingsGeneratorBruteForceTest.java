package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.generation;

import edu.rice.cs.bioinfo.library.programming.Container;
import junit.framework.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/9/13
 * Time: 5:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class AllBranchingsGeneratorBruteForceTest
{
    @Test
    public void testGenerate()
    {
        final Container<Integer> numBranchingsGeneratedAccum = new Container<Integer>(0);
        final Set<Set<String>> expectedBranchings = new HashSet<Set<String>>();
        expectedBranchings.add(new HashSet<String>(Arrays.asList("AC", "CD")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("AC", "ED")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("AC", "FD")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("BC", "CD")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("BC", "ED")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("BC", "FD")));

        new AllBranchingsGeneratorBruteForce<String>(new HashSet<String>(Arrays.asList("AC", "BC", "CD", "ED", "FD")))
        {
            @Override
            protected void branchingGenerated(Set<String> edges) {
                if(expectedBranchings.contains(edges))
                {
                    expectedBranchings.remove(edges);
                }
                else
                {
                    throw new RuntimeException("Unexpected branching");
                }
                numBranchingsGeneratedAccum.setContents(numBranchingsGeneratedAccum.getContents() + 1);
            }

            @Override
            protected Object getDestination(String edge) {
                return edge.charAt(1);
            }

            @Override
            protected Object getSource(String edge) {
                return edge.charAt(0);  //To change body of implemented methods use File | Settings | File Templates.
            }
        }.generate();

        Assert.assertEquals(new Integer(6), numBranchingsGeneratedAccum.getContents());
    }

}
