package edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.unittests;

import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.JungGraphToGraphAdapterBase;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAdditionTest;
import edu.rice.cs.bioinfo.library.programming.Func1Identity;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/9/12
 * Time: 1:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReticulateEdgeAdditionTestJung extends ReticulateEdgeAdditionTest<JungGraphToGraphAdapterBase<String, Tuple<String, String>>>
{
    @Override
    protected JungGraphToGraphAdapterBase<String, Tuple<String, String>> makeNetwork(String... nodes) {
        DirectedGraph<String, Tuple<String,String>> g = new DirectedSparseGraph<String, Tuple<String, String>>();

        for(String node : nodes)
        {
            g.addVertex(node);
        }

        return new DirectedGraphToGraphAdapter<String, Tuple<String, String>>(g, new Func1Identity<Tuple<String, String>>());
    }

    @Override
    protected boolean containsEdge(JungGraphToGraphAdapterBase<String, Tuple<String, String>> network, String source, String destination)
    {
        return network.Graph.containsEdge(new Tuple(source, destination));
    }

    @Override
    protected String makeNode(JungGraphToGraphAdapterBase<String, Tuple<String, String>> network, String node) {
       return node;
    }

    @Override
    protected Tuple<String, String> makeEdge(JungGraphToGraphAdapterBase<String, Tuple<String, String>> network, String source, String destination) {
        return new Tuple<String, String>(source, destination);
    }
}
