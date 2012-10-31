package edu.rice.cs.bioinfo.library.phylogenetics.optimization.branchlength.Brent;

import edu.rice.cs.bioinfo.library.phylogenetics.optimization.branchlength.BranchLengthOptimizer;
import edu.rice.cs.bioinfo.library.programming.Func3;
import org.junit.Test;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/20/12
 * Time: 2:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class BranchLengthOptimizerBrentAlgoApacheTest {

    @Test
    public void testOptimizeFindMaxQuad()
    {
        Func3<Object, Object, Double, Double> computeScore = new Func3<Object, Object, Double, Double>() {
            public Double execute(Object arg1, Object arg2, Double branchLength) {

                return (1.0 / (Math.pow(branchLength, 2) - branchLength + 3.0));

            }
        };

        BranchLengthOptimizerBrentAlgoApache<Object,Object,Double,Double> optimizer =
                new BranchLengthOptimizerBrentAlgoApache<Object, Object, Double, Double>(.000000000001, .0000000000000001,
                        BranchLengthOptimizerBrentAlgoApache.DOUBLE_TO_DOUBLE,
                        BranchLengthOptimizerBrentAlgoApache.DOUBLE_TO_DOUBLE);

        BranchLengthOptimizer.BranchLengthOptimizerResult<Double,Double> result = optimizer.optimize(new Object(), new Object(), computeScore);

    }
}
