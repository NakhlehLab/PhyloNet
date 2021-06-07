package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings;

import java.util.Set;

/**
 *  A strategy for finding a maximum branching of a weighted digraph.
 *
 *  A branching is defined to be a acyclic set of edges such that each node referenced by any edge is of indegree <=1.
 *  Equivalently, a branching is a forrest of zero or more arborescences (out-trees).
 *
 * @param <E> Type of directed edges in input graphs.
 * @param <W> Type of weights of edges in input graphs.
 * @see "Gibbons, Alan. "Algorithmic Graph Theory." New York: Press Syndicate of the University of Cambridge, 1985. Print."
 */
public interface MaxBranchingSolver<E,W>
{
    /**
     * Finds a maximum branching of a weighted digraph.
     *
     * That is to say a branching such that the sum of the weights of the branching edges is maximized.
     * Note that there is always a maximum branching of any digraph, even if that branching is the empty set.
     *
 //    * @param vertices the nodes of the digraph to be examined
     * @param directedEdges the directed edges of the graph to be examined
     */
    public BranchingResult<E,W> findAMaximumBranching(Set<E> directedEdges);
}
