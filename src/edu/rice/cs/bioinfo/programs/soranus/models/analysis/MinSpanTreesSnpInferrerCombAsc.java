package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.CompleteGraphFactory;
import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.FindAllMinimumSpanningTreesCombinationsAsc;
import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.GraphDisconnectedException;
import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/16/13
 * Time: 6:11 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class MinSpanTreesSnpInferrerCombAsc<S> implements MinSpanTreesSnpInferrer<S,SequencingEdge<S>,GraphDisconnectedException>
{
    public Iterable<Set<SequencingEdge<S>>> inferMinTrees(Set<S> sequencings) throws GraphDisconnectedException
    {
        final Set<SequencingEdge<S>> completeSeqGraphEdges = new CompleteGraphFactory<S,SequencingEdge<S>>()
        {
            @Override
            public SequencingEdge<S> makeEdge(S node1, S node2)
            {
                return new SequencingEdge<S>(node1, node2, getSnpDistance(node1, node2));
            }

        }.makeCompleteGraph(sequencings);

        Iterable<Set<SequencingEdge<S>>> msts =
                new FindAllMinimumSpanningTreesCombinationsAsc<Set<S>,S,SequencingEdge<S>, Long>(sequencings)
        {

                    {
                        setup();
                    }

                    @Override
                    protected Set<S> getNodesOfGraph(Set<S> graph)
                    {
                        return graph;  //To change body of implemented methods use File | Settings | File Templates.
                    }

                    @Override
                    protected Set<SequencingEdge<S>> getEdgesOfGraph(Set<S> graph)
                    {
                        return completeSeqGraphEdges;
                    }

                    @Override
                    protected Iterable<? extends SequencingEdge<S>> getIncidentEdges(final S node, Set<S> graph)
                    {
                        return IterableHelp.filter(completeSeqGraphEdges, new Predicate1<SequencingEdge<S>>()
                        {
                            public boolean execute(SequencingEdge<S> edge)
                            {
                                return edge.Sequencing1.equals(node) || edge.Sequencing2.equals(node);
                            }
                        });
                    }

                    @Override
                    protected Tuple<? extends S, ? extends S> getNodesOfEdge(SequencingEdge<S> edge, Set<S> graph)
                    {
                        return new Tuple<S, S>(edge.Sequencing1, edge.Sequencing2);
                    }

                    @Override
                    protected Long add(Long term1, Long term2)
                    {
                        return term1 + term2;
                    }

                    @Override
                    protected Long makeZero()
                    {
                        return 0l;
                    }

                    @Override
                    protected Long getWeight(SequencingEdge<S> edge)
                    {
                        return getSnpDistance(edge.Sequencing1, edge.Sequencing2);
                    }
                };

        return msts;
    }

    protected abstract Long getSnpDistance(S sequencing1, S sequencing2);
}
