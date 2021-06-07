package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.programs.soranus.viewModels.TransMapVM;
import y.base.Node;
import y.layout.tree.TreeLayouter;
import y.view.*;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/6/13
 * Time: 4:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class TransMapViewJFiles<N,E> extends TransMapViewBase
{
    public TransMapViewJFiles(TransMapVM<N,E> transMapVM) {
        super(transMapVM);
        this.setLayout(new BorderLayout());
        setup();


    }

    private void setup()
    {

        Graph2DView view = new Graph2DView();
        view.setFitContentOnResize(true);
        this.add( new ZoomableGraph2DView(view));

        Graph2D displayGraph = view.getGraph2D();
        displayGraph.getDefaultEdgeRealizer().setArrow(Arrow.STANDARD);


        TransMapVM<N,E> transMapVM = this.getTransMapVM();
        Map<N,Node> vmNodeToDisplayGraphNode = new HashMap<N, Node>();

        ShapeNodeRealizer snr = new ShapeNodeRealizer(ShapeNodeRealizer.ROUND_RECT);
        snr.setSize(100, 50);
        snr.setFillColor(new Color(245, 152, 157));


        displayGraph.setDefaultNodeRealizer(snr);


        for(N node : transMapVM.getNodes())
        {
            Node displayGraphNode = displayGraph.createNode();
            NodeLabel nodeLabel = displayGraph.getRealizer(displayGraphNode).getLabel();
            nodeLabel.setText(transMapVM.getNodeLabel(node));
            nodeLabel.setFontStyle(Font.BOLD);
            vmNodeToDisplayGraphNode.put(node, displayGraphNode);
        }

        for(E edge : transMapVM.getEdges())
        {
            N source = transMapVM.getSource(edge);
            N dest = transMapVM.getDestination(edge);

            displayGraph.createEdge(vmNodeToDisplayGraphNode.get(source), vmNodeToDisplayGraphNode.get(dest));
        }


        view.applyLayout(new TreeLayouter());

        view.fitContent();
        view.updateView();

    }
}
