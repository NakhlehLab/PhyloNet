package edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.unittests;

import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.JungGraphToGraphAdapterBase;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.UndirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree.nni.NearestNeighborInterchangeTest;
import edu.rice.cs.bioinfo.library.programming.Func1Identity;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

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

        return new UndirectedGraphToGraphAdapter<String, Tuple<String, String>>(g, new Func1Identity<Tuple<String, String>>());
    }

    @Override
    protected JungGraphToGraphAdapterBase<String, Tuple<String, String>> makeRootedTree(String... nodes)
    {
        DirectedGraph<String, Tuple<String,String>> g = new DirectedSparseGraph<String, Tuple<String, String>>();

        for(String node : nodes)
        {
            g.addVertex(node);
        }

        return new DirectedGraphToGraphAdapter<String, Tuple<String, String>>(g, new Func1Identity<Tuple<String, String>>());
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
