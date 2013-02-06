package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.programs.soranus.viewModels.TransMapVM;
import edu.uci.ics.jung.algorithms.layout.CircleLayout;
import edu.uci.ics.jung.algorithms.layout.Layout;
import edu.uci.ics.jung.algorithms.layout.ISOMLayout;
import edu.uci.ics.jung.algorithms.layout.TreeLayout;
import edu.uci.ics.jung.graph.DelegateTree;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.visualization.BasicVisualizationServer;
import org.apache.commons.collections15.Transformer;

import javax.swing.*;
import java.awt.*;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/4/13
 * Time: 7:00 PM
 * To change this template use File | Settings | File Templates.
 */
public class TransMapView<N,E> extends JPanel
{
    private final TransMapVM<N,E> _transMapVM;

    public TransMapView(TransMapVM<N, E> transMapVM) {
        _transMapVM = transMapVM;

        DelegateTree<N,E> transMap = new DelegateTree<N, E>();


        Set<N> potentialRoots = new HashSet<N>(_transMapVM.getNodes());
        for(E edge : _transMapVM.getEdges())
        {
        //    N parent = _transMapVM.getSource(edge);
            N child =  _transMapVM.getDestination(edge);
            potentialRoots.remove(child);
        }

        transMap.setRoot(potentialRoots.iterator().next());

        Set<E> edgesToAdd = new HashSet<E>(_transMapVM.getEdges());

        while(!edgesToAdd.isEmpty())
        {
            Iterator<E> edgesToAddElements = edgesToAdd.iterator();

            while(edgesToAddElements.hasNext())
            {
                E edge = edgesToAddElements.next();
                N parent = _transMapVM.getSource(edge);
                N child =  _transMapVM.getDestination(edge);
                if(transMap.getVertices().contains(parent))
                {
                    transMap.addEdge(edge, parent, child);
                    edgesToAddElements.remove();
                }
            }
        }


        Layout<N, E> layout = new TreeLayout<N,E>(transMap);
 //       layout.setSize(new Dimension(300,300));
        BasicVisualizationServer<N,E> vv = new BasicVisualizationServer<N,E>(layout);
        vv.getRenderContext().setVertexLabelTransformer(new Transformer<N,String>()
        {
            public String transform(N n) {
                return _transMapVM.getNodeLabel(n);
            }
        });
        this.add(vv);



    }
}
