package edu.rice.cs.bioinfo.library.phylogenetics.optimization.branchlength;

import edu.rice.cs.bioinfo.library.programming.Func3;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/20/12
 * Time: 1:44 PM
 * To change this template use File | Settings | File Templates.
 */
public interface BranchLengthOptimizer<G, E, L, S>
{
    public static class BranchLengthOptimizerResult<L,S>
    {
        public final L BranchLength;

        public final S Score;

        public BranchLengthOptimizerResult(L branchLength, S score)
        {
            BranchLength = branchLength;
            Score = score;
        }
    }

    public BranchLengthOptimizerResult<L,S> optimize(G graph, E edge, Func3<G,E,L,S> computeScore);
}
