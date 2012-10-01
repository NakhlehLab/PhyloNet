package edu.rice.cs.bioinfo.library.phylogenetics;

import com.sun.org.apache.bcel.internal.generic.NEW;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import junit.framework.Assert;
import org.junit.Test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
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
    class GraphMap implements GraphReadOnly<String, Tuple<String,String>>
    {
        private Map<String,HashSet<String>> _nodeToNeighbors;

        private boolean _isRooted;

        public GraphMap(Map<String,HashSet<String>> nodeToNeighbors, boolean isRooted)
        {
            _nodeToNeighbors = nodeToNeighbors;
            _isRooted = isRooted;
        }

        @Override
        public Iterable<String> getNodes() {
            return _nodeToNeighbors.keySet();
        }

        @Override
        public Iterable<Tuple<String, String>> getEdges() {

            LinkedList<Tuple<String,String>> edges = new LinkedList<Tuple<String, String>>();

            for(Map.Entry<String,HashSet<String>> kvp : _nodeToNeighbors.entrySet())
            {
                for(String dest : kvp.getValue())
                {
                    edges.add(new Tuple<String, String>(kvp.getKey(), dest));
                }
            }

            return edges;
        }

        @Override
        public Tuple<String, String> getNodesOfEdge(Tuple<String, String> edge) {
            return edge;
        }

        @Override
        public Iterable<Tuple<String, String>> getIncidentEdges(String node)
        {
            LinkedList<Tuple<String,String>> edges = new LinkedList<Tuple<String, String>>();

            for(Map.Entry<String,HashSet<String>> kvp : _nodeToNeighbors.entrySet())
            {
                for(String dest : kvp.getValue())
                {
                    if(kvp.getKey().equals(node) || dest.equals(node))
                    {
                        edges.add(new Tuple<String, String>(kvp.getKey(), dest));
                    }
                }
            }

            return edges;
        }

        @Override
        public Tuple<String, String> getEdge(String source, String destination)
        {
            return new Tuple<String, String>(source, destination);
        }

        @Override
        public boolean isRooted() {
            return _isRooted;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }

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

        try
        {
            new GraphValidator().assertValidGraph(new GraphMap(diGraph, true));
            Assert.fail("Expected exception.");
        }
        catch(IllegalArgumentException e)
        {

        }

    }


}
