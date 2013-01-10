package edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple;

import edu.rice.cs.bioinfo.library.programming.Container;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
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
public class AllDAGsGeneratorTest
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



        AllDAGsGenerator generator = new AllDAGsGenerator<String>(nodes)
        {
            @Override
            protected void graphGenerated(Iterable<Tuple<String, String>> edgesOfGraph) {

                List<Tuple<String, String>> edgesOfGraphList = IterableHelp.toList(edgesOfGraph);

                Iterator<LinkedList<Tuple<String,String>>> expectedGraphsElements = expectedGraphs.iterator();

                boolean removed = false;
                while(expectedGraphsElements.hasNext())
                {
                    LinkedList<Tuple<String,String>> expectedGraph = expectedGraphsElements.next();

                    if(edgesOfGraphList.containsAll(expectedGraph) && expectedGraph.containsAll(edgesOfGraphList))
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
        };

        generator.generate();
        Assert.assertEquals(0, expectedGraphs.size());

    }

    @Test
    public void testGenerate2()
    {
        HashSet<String> nodes = new HashSet<String>(Arrays.asList("A", "B", "C"));
        final Container<Integer> numGenerated = new Container<Integer>(0);
        AllDAGsGenerator generator = new AllDAGsGenerator<String>(nodes)
        {
            @Override
            protected void graphGenerated(Iterable<Tuple<String, String>> edgesOfGraph) {
                numGenerated.setContents(numGenerated.getContents() + 1);
            }
        };
        generator.generate();
        Assert.assertEquals(new Integer(64), numGenerated.getContents());
    }
}
