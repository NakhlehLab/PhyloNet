package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.NeighborJoiningVM;
import edu.uci.ics.jung.algorithms.layout.*;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.visualization.BasicVisualizationServer;
import org.apache.commons.collections15.Transformer;

import javax.swing.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/5/13
 * Time: 6:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class NeighborJoiningView<N,E> extends JPanel
{
    private NeighborJoiningVM<N,E> _vm;

    public NeighborJoiningView(NeighborJoiningVM<N,E> vm) {
        _vm = vm;

        UndirectedSparseGraph<N,E> graph = new UndirectedSparseGraph<N, E>();

        for(E edge : vm.Edges)
        {
            Tuple<N,N> nodes = _vm.getNodesOfEdge(edge);
            graph.addEdge(edge, nodes.Item1, nodes.Item2);
        }

        Layout<N, E> layout = new KKLayout<N, E>(graph);
        BasicVisualizationServer<N,E> vv = new BasicVisualizationServer<N,E>(layout);
        vv.getRenderContext().setVertexLabelTransformer(new Transformer<N,String>()
        {
            public String transform(N n) {
                return _vm.getNodeLabel(n);
            }
        });
        this.add(vv);
    }
}
