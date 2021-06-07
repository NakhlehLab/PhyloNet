package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.FindAllPredecessors;
import edu.rice.cs.bioinfo.library.phylogenetics.FindSuccessors;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.NodeInjector;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.NetworkValidatorBase;
import edu.rice.cs.bioinfo.library.programming.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/7/12
 * Time: 6:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReticulateEdgeAdditionInPlace<G extends Graph<N,E>,N,E> extends NetworkValidatorBase<N,E> implements ReticulateEdgeAddition<G,N,E>
{
    private final Func3<G, N, N, E> _makeEdge;

    private final Func1<G, N> _makeNode;

    private Map<N, HashSet<N>> _nodeToAncestors = new HashMap<N, HashSet<N>>();

    private Func1<G, Map<N,Set<N>>> _findAllAncestorsStrategy = new FindAllPredecessors<G,N,E>();

    public void setFindAllAncestorsStrategy(Func1<G, Map<N,Set<N>>> newStrategy)
    {
        _findAllAncestorsStrategy = newStrategy;
    }

    public ReticulateEdgeAdditionInPlace(Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge)
    {
        _makeEdge = makeEdge;
        _makeNode = makeNode;
    }

    @Override
    public Map<N,Set<N>> computeRearrangements(G network, boolean valdiateNetwork, Func4<G,E,E,E,Boolean> rearrangementComputed)
    {
        if(valdiateNetwork)
        {
            this.assertValidNetwork(network);
        }
        Map<N,Set<N>> nodeToAncestors = _findAllAncestorsStrategy.execute(network);
        computeRearrangements(network, false, rearrangementComputed, nodeToAncestors);
        return nodeToAncestors;
    }

    @Override
    public <S extends Set<N>> void computeRearrangements(G network, boolean validateNetwork, Func4<G,E,E,E,Boolean> rearrangementComputed, Map<N,S> nodeToAncestors)
    {
        if(validateNetwork)
        {
            this.assertValidNetwork(network);
        }

        N sourceEdgeGlueNode = _makeNode.execute(network);
        network.addNode(sourceEdgeGlueNode);
        N destinationEdgeGlueNode = _makeNode.execute(network);
        network.addNode(destinationEdgeGlueNode);

        LinkedList<E> networkEdges = new LinkedList<E>(); // to prevent concurrent modification exception

        for(E edge : network.getEdges())
        {
            networkEdges.add(edge);
        }

        boolean continueNeighborhoodGeneration = true;
        for(E sourceEdge : networkEdges)
        {
            if(!continueNeighborhoodGeneration)
                break;

            Tuple<N,N> nodesOfSourceEdge = network.getNodesOfEdge(sourceEdge);
            NodeInjector.NodeInjectorUndoAction<G,N,E> undoSource = NodeInjector.injectNodeIntoEdge(network, sourceEdge, sourceEdgeGlueNode, _makeEdge, false);

            for(E destinationEdge : networkEdges)
            {
                if(!continueNeighborhoodGeneration)
                    break;

                if(destinationEdge.equals(sourceEdge))
                {
                    continue;
                }

                Tuple<N,N> nodesOfDestinationEdge = network.getNodesOfEdge(destinationEdge);

                boolean wouldIntroduceCycle = wouldIntroduceCycle(nodeToAncestors, nodesOfSourceEdge, nodesOfDestinationEdge);

                if(!wouldIntroduceCycle)
                {
                    NodeInjector.NodeInjectorUndoAction<G,N,E> undoDestination = NodeInjector.injectNodeIntoEdge(network, destinationEdge, destinationEdgeGlueNode, _makeEdge, false);

                    E reticulateEdge = _makeEdge.execute(network, sourceEdgeGlueNode, destinationEdgeGlueNode);
                    network.addEdge(reticulateEdge);

                    continueNeighborhoodGeneration = rearrangementComputed.execute(network, sourceEdge, destinationEdge, reticulateEdge);

                    network.removeEdge(reticulateEdge);

                    undoDestination.undoInjection();
                }
            }

            undoSource.undoInjection();
        }
        network.removeNode(sourceEdgeGlueNode);
        network.removeNode(destinationEdgeGlueNode);
    }

    private <S extends Set<N>> boolean wouldIntroduceCycle(Map<N, S> nodeToAncestors, Tuple<N, N> nodesOfSourceEdge, Tuple<N, N> nodesOfDestinationEdge) {

        return nodesOfSourceEdge.Item1.equals(nodesOfDestinationEdge.Item2) ||
                nodeToAncestors.get(nodesOfSourceEdge.Item1).contains(nodesOfDestinationEdge.Item2);
    }

    public G performRearrangement(G network, boolean validateNetwork, E sourceEdge, E destinationEdge)
    {
        Map<N,Set<N>> nodeToAncestors = _findAllAncestorsStrategy.execute(network);
        Func<Set<N>> makeSet = new Func<Set<N>>() {
            @Override
            public Set<N> execute() {
                return new HashSet<N>();
            }
        };
        return performRearrangement(network, validateNetwork, sourceEdge, destinationEdge, nodeToAncestors, makeSet);
    }

    public <S extends Set<N>> G performRearrangement(G network, boolean validateNetwork, E sourceEdge, E destinationEdge,
                                                     Map<N,S> nodeToAncestors, Func<S> makeSet)
    {
        if(validateNetwork)
        {
            this.assertValidNetwork(network);
        }

        Tuple<N,N> nodesOfSourceEdge = network.getNodesOfEdge(sourceEdge);
        Tuple<N,N> nodesOfDestinationEdge = network.getNodesOfEdge(destinationEdge);

        boolean wouldIntroduceCycle = wouldIntroduceCycle(nodeToAncestors, nodesOfSourceEdge, nodesOfDestinationEdge);

        if(!wouldIntroduceCycle)
        {
            N sourceEdgeGlueNode = _makeNode.execute(network);
            network.addNode(sourceEdgeGlueNode);
            N destinationEdgeGlueNode = _makeNode.execute(network);
            network.addNode(destinationEdgeGlueNode);

            NodeInjector.injectNodeIntoEdge(network, sourceEdge, sourceEdgeGlueNode, _makeEdge, false);
            NodeInjector.injectNodeIntoEdge(network, destinationEdge, destinationEdgeGlueNode, _makeEdge, false);

            E reticulateEdge = _makeEdge.execute(network, sourceEdgeGlueNode, destinationEdgeGlueNode);
            network.addEdge(reticulateEdge);

            nodeToAncestors.put(sourceEdgeGlueNode, makeSet.execute());
            S destinationEdgeGlueNodeAncestors = makeSet.execute();
            nodeToAncestors.put(destinationEdgeGlueNode, destinationEdgeGlueNodeAncestors);

            updateAncestorsToIncludeNewParent(nodeToAncestors, sourceEdgeGlueNode, nodesOfSourceEdge.Item1);
            updateAncestorsToIncludeNewParent(nodeToAncestors, destinationEdgeGlueNode, sourceEdgeGlueNode);

            updateAncestorsToIncludeNewParent(nodeToAncestors, destinationEdgeGlueNode, nodesOfDestinationEdge.Item1);


            FindSuccessors<N,E> findSuccessors = new FindSuccessors<N,E>();

            for(N node : findSuccessors.execute(network, sourceEdgeGlueNode))
            {
                nodeToAncestors.get(node).add(sourceEdgeGlueNode);
            }

            for(N node : findSuccessors.execute(network, destinationEdgeGlueNode))
            {
                S ancestors = nodeToAncestors.get(node);
                ancestors.add(destinationEdgeGlueNode);
                ancestors.addAll(destinationEdgeGlueNodeAncestors);
            }

             this.assertValidNetwork(network);


            return network;
        }
        else
        {
            throw new IllegalArgumentException("Given rearrangement would introduce a cycle.");
        }

    }

    private <S extends Set<N>> void updateAncestorsToIncludeNewParent(Map<N,S> nodeToAncestors, N node, N newParent)
    {
        S parentAncestors = nodeToAncestors.get(newParent);
        S nodeAncestors = nodeToAncestors.get(node);
        nodeAncestors.addAll(parentAncestors);
        nodeAncestors.add(newParent);
    }


}
