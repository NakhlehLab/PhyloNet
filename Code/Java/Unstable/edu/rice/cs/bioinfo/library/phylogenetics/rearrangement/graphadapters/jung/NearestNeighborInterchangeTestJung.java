package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.graphadapters.jung;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree.nni.NearestNeighborInterchangeTest;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.uci.ics.jung.graph.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/6/12
 * Time: 5:24 PM
 * To change this template use File | Settings | File Templates.
 */
public class NearestNeighborInterchangeTestJung extends NearestNeighborInterchangeTest<JungGraphToGraphAdapterBase<String, Tuple<String, String>>>
{
    public NearestNeighborInterchangeTestJung()
    {
        super(new Func2<String, String, Tuple<String, String>>() {
            public Tuple<String, String> execute(String node1, String node2) {
                return new Tuple<String, String>(node1, node2);
            }
        });
    }

    @Override
    protected JungGraphToGraphAdapterBase<String, Tuple<String, String>> makeUnrootedTree(String... nodes)
    {
        Graph<String, Tuple<String,String>> g = new UndirectedSparseGraph<String, Tuple<String, String>>();

        for(String node : nodes)
        {
            g.addVertex(node);
        }

        return new UndirectedGraphToGraphAdapter<String, Tuple<String, String>>(g, new Func1<Tuple<String, String>,Tuple<String, String>>()
        {
            public Tuple<String, String> execute(Tuple<String, String> edge) {
                return edge;  //To change body of implemented methods use File | Settings | File Templates.
            }
        });
    }

    @Override
    protected JungGraphToGraphAdapterBase<String, Tuple<String, String>> makeRootedTree(String... nodes)
    {
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
    protected boolean containsEdge(JungGraphToGraphAdapterBase<String, Tuple<String, String>> tree, String source, String destination, boolean directed) {
        if(directed)
        {
            return tree.Graph.containsEdge(new Tuple<String, String>(source, destination));
        }
        else
        {
            return tree.Graph.containsEdge(new Tuple<String, String>(source, destination)) ||
                   tree.Graph.containsEdge(new Tuple<String, String>(destination, source));
        }
    }
}
