package edu.rice.cs.bioinfo.library.phylogenetics.optimization.branchlength.Brent;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/20/12
 * Time: 1:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class BranchLengthOptimizerBrentAlgoApache<G, E, L, S> implements BranchLengthOptimizerBrentAlgo<G, E, L, S>
{
    public static final Func1<Double,Double> DOUBLE_TO_DOUBLE = new Func1<Double, Double>() {
        public Double execute(Double input) {
            return input;
        }
    };

    private final Func1<S,Double> _scoreToDouble;

    private final Func1<Double,L> _doubleToLength;

    private final double _rel;

    private final double _abs;

    BranchLengthOptimizerBrentAlgoApache(double rel, double abs, Func1<Double,L> doubleToLength, Func1<S,Double> scoreToDouble)
    {
        _rel = rel;
        _abs = abs;
        _doubleToLength = doubleToLength;
        _scoreToDouble = scoreToDouble;
    }

    public BranchLengthOptimizerResult<L,S> optimize(final G graph, final E edge, final Func3<G, E, L, S> computeScore)
    {
        UnivariateFunction functionToOptimize = new UnivariateFunction() {
            public double value(double x) {
                L suggestedBranchLength = _doubleToLength.execute(x);
                S y = computeScore.execute(graph, edge, suggestedBranchLength);
                return _scoreToDouble.execute(y);
            }
        };

        BrentOptimizer optimizer = new BrentOptimizer(_rel,_abs);
        UnivariatePointValuePair maxFoundValue = optimizer.optimize(Integer.MAX_VALUE, functionToOptimize, GoalType.MAXIMIZE, 0, Double.MAX_VALUE);

        L maxFoundBL = _doubleToLength.execute(maxFoundValue.getPoint());
        S scoreForBL = computeScore.execute(graph, edge, maxFoundBL);

        return new BranchLengthOptimizerResult(maxFoundBL, scoreForBL);

    }
}
