package edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import junit.framework.Assert;
import org.junit.Test;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/9/13
 * Time: 3:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class AllDigraphsGeneratorTest
{
    @Test
    public void testGenerate()
    {
        HashSet<String> nodes = new HashSet<String>(Arrays.asList("A", "B"));

        final LinkedList<LinkedList<Tuple<String,String>>> expectedGraphs = new LinkedList<LinkedList<Tuple<String, String>>>();
        expectedGraphs.add(new LinkedList<Tuple<String, String>>());
        expectedGraphs.add(new LinkedList<Tuple<String, String>>(Arrays.asList(new Tuple<String, String>("A", "B"))));
        expectedGraphs.add(new LinkedList<Tuple<String, String>>(Arrays.asList(new Tuple<String, String>("B", "A"))));
        expectedGraphs.add(new LinkedList<Tuple<String, String>>(Arrays.asList(new Tuple<String, String>("A", "B"), new Tuple<String, String>("B", "A"))));


        AllDiGraphsGenerator<String> allDags = new AllDiGraphsGenerator<String>(nodes);

        for(Set<Tuple<String,String>> edgesOfDag : allDags)
        {
            Iterator<LinkedList<Tuple<String,String>>> expectedGraphsElements = expectedGraphs.iterator();

            boolean removed = false;
            while(expectedGraphsElements.hasNext())
            {
                LinkedList<Tuple<String,String>> expectedGraph = expectedGraphsElements.next();

                if(edgesOfDag.containsAll(expectedGraph) && expectedGraph.containsAll(edgesOfDag))
                {
                    expectedGraphsElements.remove();
                    removed = true;
                    break;
                }
            }

            if(!removed)
            {
                throw new RuntimeException("No match was found.");
            }
        }

        Assert.assertEquals(0, expectedGraphs.size());

    }

    @Test
    public void testGenerate2()
    {
        HashSet<String> nodes = new HashSet<String>(Arrays.asList("A", "B", "C"));
        int numGenerated = 0;
        AllDiGraphsGenerator<String> allDags = new AllDiGraphsGenerator<String>(nodes);

        for(Set<Tuple<String,String>> edgesOfDag : allDags)
        {
            numGenerated++;
        }
        Assert.assertEquals(64, numGenerated);
    }
}
