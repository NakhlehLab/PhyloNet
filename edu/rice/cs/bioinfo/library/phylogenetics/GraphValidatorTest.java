package edu.rice.cs.bioinfo.library.phylogenetics;

import junit.framework.Assert;
import org.junit.Test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/18/12
 * Time: 4:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphValidatorTest
{

    @Test
    public void testAssertValidGraphWithCycle1()
    {
        Map<String, HashSet<String>> diGraph = new HashMap<String, HashSet<String>>();

        HashSet<String> aDest = new HashSet<String>();
        aDest.add("B");

        HashSet<String> bDest = new HashSet<String>();
        bDest.add("C");

        HashSet<String> cDest = new HashSet<String>();
        cDest.add("A");

        diGraph.put("A", aDest);
        diGraph.put("B", bDest);
        diGraph.put("C", cDest);

        PhyloGraph<String> graph = new DirectedPhyloGraphDefault<String>();
        graph.addNodes("A", "B", "C");
        graph.addEdges(new PhyloEdge<String>("A", "B"), new PhyloEdge<String>("B", "C"), new PhyloEdge<String>("C", "A"));

        try
        {
            new GraphValidator().assertValidGraph(graph);
            Assert.fail("Expected exception.");
        }
        catch(IllegalArgumentException e)
        {

        }

    }


}
