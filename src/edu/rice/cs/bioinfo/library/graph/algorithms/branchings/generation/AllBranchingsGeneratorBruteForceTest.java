package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.generation;

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
       // final Container<Integer> numBranchingsGeneratedAccum = new Container<Integer>(0);
        final Set<Set<String>> expectedBranchings = new HashSet<Set<String>>();
        expectedBranchings.add(new HashSet<String>());
        expectedBranchings.add(new HashSet<String>(Arrays.asList("ED")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("FD")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("CD")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("AC")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("BC")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("AC", "CD")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("AC", "ED")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("AC", "FD")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("BC", "CD")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("BC", "ED")));
        expectedBranchings.add(new HashSet<String>(Arrays.asList("BC", "FD")));

        AllBranchingsGeneratorBruteForce<String> allBranchings = new AllBranchingsGeneratorBruteForce<String>(new HashSet<String>(Arrays.asList("AC", "BC", "CD", "ED", "FD")))
        {
            @Override
            protected Object getSource(String edge) {
                return edge.charAt(0);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Object getDestination(String edge) {
                return edge.charAt(1);
            }
        };


        for(Set<String> branching : allBranchings)
        {
            if(expectedBranchings.contains(branching))
            {
                expectedBranchings.remove(branching);
            }
            else
            {
                throw new RuntimeException("Unexpected branching");
            }
        }

        Assert.assertEquals(0, expectedBranchings.size());
    }

}
