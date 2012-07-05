package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/27/12
 * Time: 2:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class AddReticulationEdge<N,E> implements Proc6<Graph<N,E>, E, E, Func1<Graph<N,E>,N>, Func3<Graph<N,E>,N,N,E>, Boolean>
{
    public void execute(Graph<N,E> graph, E edge1, E edge2, Func1<Graph<N,E>,N> makeNode, Func3<Graph<N,E>,N,N,E> makeEdge, Boolean performValidation)
    {
        N node1 = makeNode.execute(graph);
        graph.addNode(node1);
        NodeInjector.NodeInjectorUndoAction<Graph<N,E>,N,E> undo1 = NodeInjector.injectNodeIntoEdge(graph, edge1, node1, makeEdge, false);

        N node2 = makeNode.execute(graph);
        graph.addNode(node2);
        NodeInjector.NodeInjectorUndoAction<Graph<N,E>,N,E> undo2 = NodeInjector.injectNodeIntoEdge(graph, edge2, node2, makeEdge, false);

        E reticulateEdge = makeEdge.execute(graph, node1, node2);
        graph.addEdge(reticulateEdge);

        if(performValidation)
        {
            try
            {
                GraphValidator.assertValidGraph(graph);
            }
            catch(IllegalArgumentException e)
            {
                graph.removeEdge(reticulateEdge);
                undo2.undoInjection();
                graph.removeNode(node2);
                undo1.undoInjection();
                graph.removeNode(node1);
                throw e;
            }
        }

    }
}
