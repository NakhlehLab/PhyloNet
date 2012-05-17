package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.*;
import edu.rice.cs.bioinfo.library.programming.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/7/12
 * Time: 6:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReticulateEdgeAdditionInPlace<G extends Graph<N,E>,N,E> extends ReticulateEdgeAdditionBase<G,N,E>
{
    private final Func3<G, N, N, E> _makeEdge;

    private final Func1<G, N> _makeNode;

    public ReticulateEdgeAdditionInPlace(Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge)
    {
        _makeEdge = makeEdge;
        _makeNode = makeNode;
    }

    public void computeRearrangementsWithoutValidation(G network, Proc4<G,E,E,E> rearrangementComputed)
    {
        N sourceEdgeGlueNode = _makeNode.execute(network);
        network.addNode(sourceEdgeGlueNode);
        N destinationEdgeGlueNode = _makeNode.execute(network);
        network.addNode(destinationEdgeGlueNode);

        LinkedList<E> networkEdges = new LinkedList<E>(); // to prevent concurrent modification exception

        for(E edge : network.getEdges())
        {
            networkEdges.add(edge);
        }

        for(E sourceEdge : networkEdges)
        {
            Tuple<N,N> nodesOfSourceEdge = network.getNodesOfEdge(sourceEdge);
            NodeInjector.NodeInjectorUndoAction<G,N,E> undoSource = NodeInjector.injectNodeIntoEdge(network, sourceEdge, sourceEdgeGlueNode, _makeEdge, false);

            for(E destinationEdge : networkEdges)
            {
                if(destinationEdge.equals(sourceEdge))
                {
                    continue;
                }

                Tuple<N,N> nodesOfDestinationEdge = network.getNodesOfEdge(destinationEdge);



                if(!isPath(network, nodesOfDestinationEdge.Item2, nodesOfSourceEdge.Item1))
                {
                    NodeInjector.NodeInjectorUndoAction<G,N,E> undoDestination = NodeInjector.injectNodeIntoEdge(network, destinationEdge, destinationEdgeGlueNode, _makeEdge, false);

                    E reticulateEdge = _makeEdge.execute(network, sourceEdgeGlueNode, destinationEdgeGlueNode);
                    network.addEdge(reticulateEdge);

                    rearrangementComputed.execute(network, sourceEdge, destinationEdge, reticulateEdge);

                    network.removeEdge(reticulateEdge);

                    undoDestination.undoInjection();
                }
            }

            undoSource.undoInjection();
        }
        network.removeNode(sourceEdgeGlueNode);
        network.removeNode(destinationEdgeGlueNode);
    }

    public G perormRearrangementWithoutValidation(G network, E sourceEdge, E destinationEdge)
    {
        Tuple<N,N> nodesOfSourceEdge = network.getNodesOfEdge(sourceEdge);
        Tuple<N,N> nodesOfDestinationEdge = network.getNodesOfEdge(destinationEdge);

        N sourceEdgeGlueNode = _makeNode.execute(network);
        network.addNode(sourceEdgeGlueNode);
        N destinationEdgeGlueNode = _makeNode.execute(network);
        network.addNode(destinationEdgeGlueNode);

        NodeInjector.injectNodeIntoEdge(network, sourceEdge, sourceEdgeGlueNode, _makeEdge, false);

        if(!isPath(network, nodesOfDestinationEdge.Item2, nodesOfSourceEdge.Item1))
        {
            NodeInjector.injectNodeIntoEdge(network, destinationEdge, destinationEdgeGlueNode, _makeEdge, false);

            E reticulateEdge = _makeEdge.execute(network, sourceEdgeGlueNode, destinationEdgeGlueNode);
            network.addEdge(reticulateEdge);

            return network;
        }
        else
        {
            throw new IllegalArgumentException("Given rearrangement would introduce a cycle.");
        }

    }

    protected boolean isPath(G network, N start, N end)
    {
        LinkedList<N> toExpand = new LinkedList<N>();

        toExpand.push(start);

        while(toExpand.size() > 0)
        {
            N node = toExpand.pop();

            if(end.equals(node))
            {
                return true;
            }

            for(E edge : network.getIncidentEdges(node))
            {
                Tuple<N,N> nodesOfEdge = network.getNodesOfEdge(edge);

                if(nodesOfEdge.Item1.equals(node))
                {
                    toExpand.push(nodesOfEdge.Item2);
                }
            }
        }


        return false;

    }
}
