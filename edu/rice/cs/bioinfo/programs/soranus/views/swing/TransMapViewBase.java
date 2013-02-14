package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.programs.soranus.viewModels.TransMapVM;

import javax.swing.*;
import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/4/13
 * Time: 7:00 PM
 * To change this template use File | Settings | File Templates.
 */
public class TransMapViewBase<N,E> extends JPanel
{
    private final TransMapVM<N,E> _transMapVM;

    protected TransMapVM<N,E> getTransMapVM()
    {
        return _transMapVM;
    }

    public TransMapViewBase(TransMapVM<N, E> transMapVM) {
        _transMapVM = transMapVM;
        /*
        DelegateTree<N,E> transMap = new DelegateTree<N, E>();




        transMap.setRoot(findRoot());

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
        this.add(vv); */
    }

    protected N findRoot()
    {
        Set<N> potentialRoots = new HashSet<N>(_transMapVM.getNodes());
        for(E edge : _transMapVM.getEdges())
        {
            N child =  _transMapVM.getDestination(edge);
            potentialRoots.remove(child);
        }

        if(potentialRoots.size() == 1)
        {
            return potentialRoots.iterator().next();
        }
        else
        {
            throw new IllegalStateException();
        }
    }
}
