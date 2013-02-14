package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.NeighborJoiningVM;
import y.base.Edge;
import y.base.Node;
import y.layout.tree.BalloonLayouter;
import y.view.*;

import javax.swing.*;
import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/5/13
 * Time: 6:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class NeighborJoiningViewJFiles<N,E> extends NeighborJoiningViewJFilesBase<N,E>
{
    JScrollPane _scrollPane = new JScrollPane();

    public NeighborJoiningViewJFiles(NeighborJoiningVM<N, E> vm) {
        super(vm);
        this.setLayout(new BorderLayout());
        setup();




    }

    private void setup()
    {
        NeighborJoiningVM<N,E> viewModel = this.getViewModel();
        final Graph2DView view = new Graph2DView();
        view.setFitContentOnResize(true);
        this.add( new ZoomableGraph2DView(view));

        Graph2D displayGraph = view.getGraph2D();
        displayGraph.getDefaultEdgeRealizer().setArrow(Arrow.NONE);

        Map<N,Node> vmNodeToDisplayGraphNode = new HashMap<N, Node>();

        ShapeNodeRealizer leafNodeRealizer = new ShapeNodeRealizer(ShapeNodeRealizer.ROUND_RECT);
        leafNodeRealizer.setFillColor(new Color(245, 152, 157));
        ShapeNodeRealizer internalNodeRealizer = new ShapeNodeRealizer(ShapeNodeRealizer.ELLIPSE);
        internalNodeRealizer.setFillColor(new Color(122, 204,200));



        for(N node : viewModel.Nodes)
        {
            Node displayGraphNode = displayGraph.createNode();
            vmNodeToDisplayGraphNode.put(node, displayGraphNode);

            NodeRealizer nodeRealizer = viewModel.isLeaf(node) ? leafNodeRealizer.createCopy() : internalNodeRealizer.createCopy();

            displayGraph.setRealizer(displayGraphNode, nodeRealizer);

            NodeLabel nodeLabel = displayGraph.getRealizer(displayGraphNode).getLabel();
            nodeLabel.setText(viewModel.getNodeLabel(node));
            nodeLabel.setFontStyle(Font.BOLD);



        }

        for(E edge : viewModel.Edges)
        {
            Tuple<N,N> nodesOfEdge = viewModel.getNodesOfEdge(edge);

            Edge dispalyGraphEdge = displayGraph.createEdge(vmNodeToDisplayGraphNode.get(nodesOfEdge.Item1), vmNodeToDisplayGraphNode.get(nodesOfEdge.Item2));
            EdgeLabel edgeLabel = displayGraph.getRealizer(dispalyGraphEdge).getLabel();

            edgeLabel.setText(viewModel.getEdgeLabel(edge));
        }

        view.applyLayout(new BalloonLayouter());
        view.fitWorldRect();
        view.updateView();
    }
}
