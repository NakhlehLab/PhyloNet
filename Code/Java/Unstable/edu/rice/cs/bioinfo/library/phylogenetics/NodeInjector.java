package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Tuple;

public class NodeInjector
{
    public static class NodeInjectorUndoAction<G extends Graph<N,E>, N, E>
    {
        E __removedEdge;private E get_removedEdge(){ return __removedEdge; }; private void set_removedEdge(E value){__removedEdge = value; };

        E __addedEdge1;private E get_addedEdge1(){ return __addedEdge1; }; private void set_addedEdge1(E value){__addedEdge1 = value; };

        E __addedEdge2;private E get_addedEdge2(){ return __addedEdge2; }; private void set_addedEdge2(E value){__addedEdge2 = value; };

        G __graph;private G get_graph(){ return __graph; }; private void set_graph(G value){__graph = value; };

        NodeInjectorUndoAction(E removedEdge, E addedEdge1, E addedEdge2, G graph)
        {
            set_removedEdge(removedEdge);
            set_addedEdge1(addedEdge1);
            set_addedEdge2(addedEdge2);
            set_graph(graph);
        }

        public void undoInjection()
        {
            get_graph().removeEdge(get_addedEdge1());
            get_graph().removeEdge(get_addedEdge2());
            get_graph().addEdge(get_removedEdge());
        }


    }

    public static <G extends Graph<N,E>, N, E> NodeInjectorUndoAction<G,N,E> injectNodeIntoEdge(G graph, E edge, N node, Func3<G, N, N, E> makeEdge, boolean makeNodeOnlySource)
    {

        Tuple<N, N> edgeNodes = graph.getNodesOfEdge(edge);
        boolean item1isSource = graph.isRooted() &&  new IsDestinationNode().execute(graph, edgeNodes.Item2, edge);
        graph.removeEdge(edge);
        E edgeToAdd1 = null;
        E edgeToAdd2 = null;

        if(graph.isRooted() && makeNodeOnlySource)
        {
            edgeToAdd1 = makeEdge.execute(graph, node, edgeNodes.Item1);
            edgeToAdd2 = makeEdge.execute(graph, node, edgeNodes.Item2);
        }
        else if(item1isSource)
        {
            edgeToAdd1 = makeEdge.execute(graph, edgeNodes.Item1, node);
            edgeToAdd2 = makeEdge.execute(graph, node, edgeNodes.Item2);
        }
        else
        {
            edgeToAdd1 = makeEdge.execute(graph, edgeNodes.Item2, node);
            edgeToAdd2 = makeEdge.execute(graph, node, edgeNodes.Item1);
        }

        graph.addEdge(edgeToAdd1);
        graph.addEdge(edgeToAdd2);

        return new NodeInjectorUndoAction<G, N, E>(edge, edgeToAdd1, edgeToAdd2, graph);
    }
}

