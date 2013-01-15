package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/14/13
 * Time: 12:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class BranchingResult<E,W>
{
    public final Set<E> Branching;

    public final W BranchingWeight;

    public BranchingResult(Set<E> branching, W branchingWeight) {
        Branching = branching;
        BranchingWeight = branchingWeight;
    }
}