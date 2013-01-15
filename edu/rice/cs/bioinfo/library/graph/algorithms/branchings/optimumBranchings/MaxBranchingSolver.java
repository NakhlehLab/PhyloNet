package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/14/13
 * Time: 12:45 PM
 * To change this template use File | Settings | File Templates.
 */
public interface MaxBranchingSolver<V,E,W>
{
    public BranchingResult<E,W> findAMaximumBranching(Set<V> vertices, Set<E> directedEdges);
}
