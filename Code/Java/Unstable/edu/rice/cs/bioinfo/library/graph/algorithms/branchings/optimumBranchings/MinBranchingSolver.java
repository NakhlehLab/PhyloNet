package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings;

import java.util.Set;

/**
 *  A strategy for finding a minimum branching of a weighted digraph.
 *
 *  A branching is defined to be a acyclic set of edges such that each node referenced by any edge is of indegree <=1.
 *  Equivalently, a branching is a forrest of zero or more arborescences (out-trees).
 *
 * @param <E> type of directed edges in input graphs.
 * @param <W> type of weights of edges in input graphs.
 * @see "Gibbons, Alan. "Algorithmic Graph Theory." New York: Press Syndicate of the University of Cambridge, 1985. Print."
 */

public interface MinBranchingSolver<E,W>
{
    /**
     * Finds a minimum branching of a weighted digraph.
     *
     * That is to say a branching such that the sum of the weights of the branching edges is minimal.
     * Note that there is always a minimum branching of any digraph, even if that branching is the empty set.
     *
     * @param directedEdges the directed edges of the graph to be examined
     */
    public BranchingResult<E,W> findMinimumBranching(Set<E> directedEdges);

}
