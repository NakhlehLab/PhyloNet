package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.CompleteGraphFactory;
import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.GraphDisconnectedException;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/22/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class MinSpanTreesSnpInferrerMaxTwo<S>
{
    public Iterable<Set<SequencingEdge<S>>> inferMinTrees(Set<S> sequencings)
    {
        final Map<S,Set<S>> nodeToConnectedSet = makeConnectedSets(sequencings);

        Set<SequencingEdge<S>> completeSeqGraphEdges = new CompleteGraphFactory<S,SequencingEdge<S>>()
        {
            @Override
            public SequencingEdge<S> makeEdge(S node1, S node2)
            {
                return new SequencingEdge<S>(node1, node2, getSnpDistance(node1, node2));
            }

        }.makeCompleteGraph(sequencings);

        List<SequencingEdge<S>> sequencingEdgesByDistanceAsc = new ArrayList<SequencingEdge<S>>(completeSeqGraphEdges);
        Collections.sort(sequencingEdgesByDistanceAsc);
        sequencingEdgesByDistanceAsc = new LinkedList<SequencingEdge<S>>(sequencingEdgesByDistanceAsc);

        Set<SequencingEdge<S>> mst = new HashSet<SequencingEdge<S>>();
        for(int i = 0; i<sequencingEdgesByDistanceAsc.size(); i++)
        {
            SequencingEdge<S> ithEdge = sequencingEdgesByDistanceAsc.get(i);
            if(nodeToConnectedSet.get(ithEdge.Sequencing1) == nodeToConnectedSet.get(ithEdge.Sequencing2))
                continue;

            Collection<SequencingEdge<S>> preSelectConsideration =
                    findConsiderationEdges(sequencingEdgesByDistanceAsc, i, ithEdge.SnpDistance, nodeToConnectedSet);
            mst.add(ithEdge);
            unionConnectedSets(ithEdge, nodeToConnectedSet);
            Collection<SequencingEdge<S>> postSelectConsideration =
                    findConsiderationEdges(sequencingEdgesByDistanceAsc, i, ithEdge.SnpDistance, nodeToConnectedSet);

            if(preSelectConsideration.size() != postSelectConsideration.size() + 1)
            {
                preSelectConsideration.remove(ithEdge);
                preSelectConsideration.removeAll(postSelectConsideration);

                SequencingEdge<S> nullifiedOption = preSelectConsideration.iterator().next();

                return Arrays.asList(findAMstWithEdge(nullifiedOption, completeSeqGraphEdges, sequencings),
                                     findAMstWithEdge(ithEdge, completeSeqGraphEdges, sequencings));
            }

        }

        if(mst.size() != (sequencings.size() - 1))
        {
            throw new RuntimeException(new GraphDisconnectedException());
        }

        return new HashSet<Set<SequencingEdge<S>>>(Arrays.asList(mst));

    }

    private List<SequencingEdge<S>> findConsiderationEdges(List<SequencingEdge<S>> potentialEdges, int startIndex, long weight, Map<S,Set<S>> nodeToConnectedSet)
    {
        List<SequencingEdge<S>> considerationEdges = new LinkedList<SequencingEdge<S>>();

        for(int i = startIndex; i<potentialEdges.size(); i++)
        {
            SequencingEdge<S> ithEdge = potentialEdges.get(i);
            if(ithEdge.SnpDistance == weight)
            {
                if(nodeToConnectedSet.get(ithEdge.Sequencing1) != nodeToConnectedSet.get(ithEdge.Sequencing2))
                {
                    considerationEdges.add(ithEdge);
                }
            }

        }
        return considerationEdges;
    }

    protected abstract Long getSnpDistance(S node1, S node2);

    private Map<S, Set<S>> makeConnectedSets(Set<S> graphNodes)
    {
        Map<S,Set<S>> nodeToConnectedSet = new HashMap<S,Set<S>>();
        for(S node : graphNodes)
        {
            nodeToConnectedSet.put(node, new HashSet<S>(Arrays.asList(node)));
        }

        return nodeToConnectedSet;
    }

    private Set<SequencingEdge<S>> findAMstWithEdge(SequencingEdge<S> mstEdge, Set<SequencingEdge<S>> graphEdges, Set<S> graphNodes)
    {
        Map<S,Set<S>> nodeToConnectedSet = makeConnectedSets(graphNodes);

        Set<SequencingEdge<S>> mst = new HashSet<SequencingEdge<S>>();
        mst.add(mstEdge);
        unionConnectedSets(mstEdge, nodeToConnectedSet);

        List<SequencingEdge<S>> edgesByDistanceAsc = new ArrayList<SequencingEdge<S>>(graphEdges);
        Collections.sort(edgesByDistanceAsc);
        edgesByDistanceAsc = new LinkedList<SequencingEdge<S>>(edgesByDistanceAsc);

        while(!edgesByDistanceAsc.isEmpty())
        {
            SequencingEdge<S> edge = edgesByDistanceAsc.remove(0);
            if(nodeToConnectedSet.get(edge.Sequencing1) == nodeToConnectedSet.get(edge.Sequencing2))
                continue;

            mst.add(edge);
            unionConnectedSets(edge, nodeToConnectedSet);
        }

        return mst;
    }


    private void unionConnectedSets(SequencingEdge<S> edge, Map<S, Set<S>> nodeToConnectedSet)
    {
        S node1 = edge.Sequencing1;
        S node2 = edge.Sequencing2;

        Set<S> set1 = nodeToConnectedSet.get(node1);
        Set<S> set2 = nodeToConnectedSet.get(node2);

        set1.addAll(set2);

        for(S nodeOfSet2 : set2)
        {
            nodeToConnectedSet.put(nodeOfSet2, set1);
        }

    }


    private List<SequencingEdge<S>> getEdgesOfHeadDistance(List<SequencingEdge<S>> sequencingEdgesByDistanceAsc)
    {
        if(sequencingEdgesByDistanceAsc.isEmpty())
            throw new IllegalArgumentException();

        LinkedList<SequencingEdge<S>> edgesOfHeadWeight = new LinkedList<SequencingEdge<S>>();
        Long headDistance = sequencingEdgesByDistanceAsc.get(0).SnpDistance;

        for(int i = 0; i<sequencingEdgesByDistanceAsc.size() &&
                         sequencingEdgesByDistanceAsc.get(i).SnpDistance == headDistance; i++ )
        {
            edgesOfHeadWeight.add(sequencingEdgesByDistanceAsc.get(i));
        }

        return edgesOfHeadWeight;

    }


}
