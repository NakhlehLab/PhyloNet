package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.CompleteGraphFactory;
import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.FindAllMinimumSpanningTreesKruskalTies;
import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.GraphDisconnectedException;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 5/3/13
 * Time: 1:34 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class FindMinSpanTreesSnp<S>
{
    public Iterable<Set<SequencingEdge<S>>> inferMinTrees(final Set<S> sequencings) throws GraphDisconnectedException
    {
        final Map<S, Collection<SequencingEdge<S>>> nodeToIndcidentEdges = new HashMap<S,Collection<SequencingEdge<S>>>();
        for(S node : sequencings)
        {
            nodeToIndcidentEdges.put(node, new LinkedList<SequencingEdge<S>>());
        }

        final Set<SequencingEdge<S>> completeSeqGraphEdges = new CompleteGraphFactory<S,SequencingEdge<S>>()
                {
                    @Override
                    public SequencingEdge<S> makeEdge(S node1, S node2)
                    {
                        SequencingEdge<S> edge = new SequencingEdge<S>(node1, node2, getSnpDistance(node1, node2));
                        nodeToIndcidentEdges.get(node1).add(edge);
                        nodeToIndcidentEdges.get(node2).add(edge);
                        return edge;
                    }

                }.makeCompleteGraph(sequencings);

        return new FindAllMinimumSpanningTreesKruskalTies<Object,S,SequencingEdge<S>,Long>(sequencings)
        {

            @Override
            protected Long add(Long term1, Long term2)
            {
                return term1 + term1;
            }

            @Override
            protected Long makeZero()
            {
                return 0l;
            }

            @Override
            protected Set<S> getNodesOfGraph(Object graph)
            {
                return sequencings;
            }

            @Override
            protected Set<SequencingEdge<S>> getEdgesOfGraph(Object graph)
            {
                return completeSeqGraphEdges;
            }

            @Override
            protected Iterable<? extends SequencingEdge<S>> getIncidentEdges(S node, Object graph)
            {
                return nodeToIndcidentEdges.get(node);
            }

            @Override
            protected Tuple<? extends S, ? extends S> getNodesOfEdge(SequencingEdge<S> edge, Object graph)
            {
                return new Tuple<S, S>(edge.Sequencing1, edge.Sequencing2);
            }

            @Override
            protected Long getWeightOfEdge(SequencingEdge<S> edge)
            {
                return edge.SnpDistance;
            }
        };
    }

    protected abstract Long getSnpDistance(S sequencing1, S sequencing2);
}
