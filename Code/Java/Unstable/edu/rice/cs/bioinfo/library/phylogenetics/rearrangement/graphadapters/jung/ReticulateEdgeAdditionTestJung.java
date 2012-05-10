package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.graphadapters.jung;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.*;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.uci.ics.jung.graph.*;

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
        Graph<String, Tuple<String,String>> g = new DirectedSparseGraph<String, Tuple<String, String>>();

        for(String node : nodes)
        {
            g.addVertex(node);
        }

        return new DirectedGraphToGraphAdapter<String, Tuple<String, String>>(g, new Func1<Tuple<String, String>,Tuple<String, String>>()
        {
            public Tuple<String, String> execute(Tuple<String, String> edge) {
                return edge;  //To change body of implemented methods use File | Settings | File Templates.
            }
        });
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
