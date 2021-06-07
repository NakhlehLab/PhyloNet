package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.CompleteGraphFactory;
import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.FindAMinimumSpanningTreeKruskal;
import edu.rice.cs.bioinfo.library.math.discrete.Combinations;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 5/28/13
 * Time: 3:48 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class MinSpanTreesSnpInferrerDiversity<S>
{
    public Iterable<Set<SequencingEdge<S>>> inferMinTrees(Set<S> sequencings)  throws Exception
    {
        Set<SequencingEdge<S>> completeSeqGraphEdges = new CompleteGraphFactory<S,SequencingEdge<S>>()
             {
                 @Override
                 public SequencingEdge<S> makeEdge(S node1, S node2)
                 {
                     return new SequencingEdge<S>(node1, node2, getSnpDistance(node1, node2));
                 }

             }.makeCompleteGraph(sequencings);

                FindAMinimumSpanningTreeKruskal<Set<SequencingEdge<S>>,SequencingEdge<S>,Long> findMst = new FindAMinimumSpanningTreeKruskal<Set<SequencingEdge<S>>,SequencingEdge<S>,Long>()
                {

                    @Override
                    protected Long add(Long term1, Long term2)
                    {
                        return term1 + term2;  //To change body of implemented methods use File | Settings | File Templates.
                    }

                    @Override
                    protected Long makeZero()
                    {
                        return 0l;  //To change body of implemented methods use File | Settings | File Templates.
                    }

                    @Override
                    protected Tuple<?, ?> getNodesOfEdge(SequencingEdge<S> edge, Set<SequencingEdge<S>> graph)
                    {
                        return new Tuple<S, S>(edge.Sequencing1, edge.Sequencing2);
                    }

                    @Override
                    protected Long getWeight(SequencingEdge<S> edge, Set<SequencingEdge<S>> graph)
                    {
                        return edge.SnpDistance;
                    }

                    @Override
                    protected Set<? extends SequencingEdge<S>> getEdges(Set<SequencingEdge<S>> graph)
                    {
                        return graph;
                    }

                    @Override
                    protected Set<?> getNodes(Set<SequencingEdge<S>> graph)
                    {
                        HashSet<S> nodes = new HashSet<S>();
                        for(SequencingEdge<S> e : graph)
                        {
                            nodes.add(e.Sequencing1);
                            nodes.add(e.Sequencing2);
                        }

                        return nodes;


                    }
                };

        Tuple<Set<SequencingEdge<S>>,Long> mstResult = findMst.execute(completeSeqGraphEdges);

        Set<SequencingEdge<S>> someMst = mstResult.Item1;
        Long mstWeight = mstResult.Item2;

        Set<Set<SequencingEdge<S>>> showMsts = new HashSet<Set<SequencingEdge<S>>>();
        showMsts.add(someMst);

        for(int i = someMst.size(); i>=2; i--)
        {
            for(Iterable<SequencingEdge<S>> excludeEdges : new Combinations<SequencingEdge<S>>(someMst, i))
            {
                Set<SequencingEdge<S>> completeSeqGraphEdgesClone = new HashSet<SequencingEdge<S>>(completeSeqGraphEdges);
                completeSeqGraphEdgesClone.removeAll(IterableHelp.toList(excludeEdges));
                mstResult = findMst.execute(completeSeqGraphEdgesClone);

                if(mstResult.Item2.equals(mstWeight))
                {
                    showMsts.add(mstResult.Item1);
                }
            }
        }

        return showMsts;
    }

    protected abstract Long getSnpDistance(S node1, S node2);
}
