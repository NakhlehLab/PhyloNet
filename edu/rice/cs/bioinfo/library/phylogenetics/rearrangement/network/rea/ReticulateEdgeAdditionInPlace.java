package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.*;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.*;
import edu.rice.cs.bioinfo.library.programming.*;
import org.junit.internal.matchers.Each;

import java.nio.ReadOnlyBufferException;
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

                if(!isReachable(network, nodesOfSourceEdge.Item1, nodesOfDestinationEdge.Item2))
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

    protected boolean isReachable(G network, N sourceNode, N destinationNode)
    {
        HashSet<N> ancestorsOfSourceNode = new HashSet<N>();
        LinkedList<N> toExpand = new LinkedList<N>();
        toExpand.push(sourceNode);

        while(toExpand.size() > 0)
        {
            N node = toExpand.pop();

            for(E edge : network.getIncidentEdges(node))
            {
                Tuple<N,N> nodesOfEdge = network.getNodesOfEdge(edge);

                if(nodesOfEdge.Item2.equals(node))
                {
                   if(!ancestorsOfSourceNode.contains(nodesOfEdge.Item1))
                   {
                       ancestorsOfSourceNode.add(nodesOfEdge.Item1);
                       toExpand.push(nodesOfEdge.Item1);
                   }
                }
            }
        }

        toExpand.push(destinationNode);

        while(toExpand.size() > 0)
        {
            N node = toExpand.pop();

            if(ancestorsOfSourceNode.contains(node))
            {
                return true;
            }
            else
            {
                for(E edge : network.getIncidentEdges(node))
                {
                    Tuple<N,N> nodesOfEdge = network.getNodesOfEdge(edge);

                    if(nodesOfEdge.Item1.equals(node))
                    {
                        toExpand.push(nodesOfEdge.Item2);
                    }
                }
            }
        }

        return false;

    }
}
