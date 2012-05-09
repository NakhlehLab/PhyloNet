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
    public ReticulateEdgeAdditionTestJung()
    {
        super(new Func1<JungGraphToGraphAdapterBase<String, Tuple<String, String>>, String>()
        {           int i = 0;
                    public String execute(JungGraphToGraphAdapterBase<String, Tuple<String, String>> input1) {
                        return "_" + (++i);
                    }
                },
        new Func3<JungGraphToGraphAdapterBase<String, Tuple<String, String>>, String, String, Tuple<String, String>>()
        {
            public Tuple<String, String> execute(JungGraphToGraphAdapterBase<String, Tuple<String, String>> arg1, String node1, String node2) {
                return new Tuple<String, String>(node1, node2);
            }
        });
    }

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
}
