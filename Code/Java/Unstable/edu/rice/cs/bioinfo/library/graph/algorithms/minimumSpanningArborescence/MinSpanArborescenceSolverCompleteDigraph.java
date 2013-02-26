package edu.rice.cs.bioinfo.library.graph.algorithms.minimumSpanningArborescence;

import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.BranchingResult;
import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.EdmondsAlgoGibbons85;
import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.MaxBranchingSolver;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.LinkedList;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/14/13
 * Time: 6:41 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class MinSpanArborescenceSolverCompleteDigraph<E,W extends Comparable<W>> implements MinSpanArborescenceSolver<E,W>
{
    private final W _one;

    private final W _zero;

    public MinSpanArborescenceSolverCompleteDigraph(W zero, W one)
    {
        _zero = zero;
        _one = one;
    }

    private class EdgeAndWeight implements Comparable<EdgeAndWeight>
    {
        public final E Edge;

        public final W Weight;

        private EdgeAndWeight(E edge, W weight) {
            Edge = edge;
            Weight = weight;
        }

        public int compareTo(EdgeAndWeight o) {
            return this.Weight.compareTo(o.Weight);
        }
    }

    public MinSpanArborescence<E,W> tryFindMinSpanArborescence(Set<E> directedEdges)
    {
        LinkedList<EdgeAndWeight> weightedEdges = new LinkedList<EdgeAndWeight>();

        for(E edge : directedEdges)
        {
            weightedEdges.add(new EdgeAndWeight(edge, getWeightOfEdge(edge)));
        }

        EdgeAndWeight maxWeightedEdge = IterableHelp.maxes(weightedEdges).iterator().next();
        W maxEdgeWeight = maxWeightedEdge.Weight;

        final W perBranchQuota = addWeight(maxEdgeWeight, _one);

        MaxBranchingSolver<E,W> maxBranchSolver = new EdmondsAlgoGibbons85<E,W>(_zero)
        {

            @Override
            protected W addWeight(W w1, W w2) {
                return MinSpanArborescenceSolverCompleteDigraph.this.addWeight(w1,w2);
            }

            @Override
            protected W subtractWeight(W w1, W w2) {
                return MinSpanArborescenceSolverCompleteDigraph.this.subtractWeight(w1, w2);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected W getWeightOfEdge(E edge) {
                W originalEdgeWeight = MinSpanArborescenceSolverCompleteDigraph.this.getWeightOfEdge(edge);
                W penality = subtractWeight(perBranchQuota, originalEdgeWeight);
                return penality;
            }

            @Override
            protected Object getSource(E edge) {
                return MinSpanArborescenceSolverCompleteDigraph.this.getSource(edge);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Object getDestination(E edge) {
                return MinSpanArborescenceSolverCompleteDigraph.this.getDestination(edge);  //To change body of implemented methods use File | Settings | File Templates.
            }
        };

        BranchingResult<E, W> br = maxBranchSolver.findAMaximumBranching(directedEdges);
        W spanWeightAccum = _zero;
        for(E spanEdge : br.Branching)
        {
            spanWeightAccum = addWeight(spanWeightAccum, getWeightOfEdge(spanEdge));
        }

        return new MinSpanArborescence<E, W>(br.Branching, spanWeightAccum);

    }

    protected abstract W addWeight(W w1, W w2);

    protected abstract W subtractWeight(W w1, W w2);

    protected abstract W getWeightOfEdge(E edge);

    protected abstract Object getSource(E edge);

    protected abstract Object getDestination(E edge);

}
