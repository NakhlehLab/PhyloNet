package edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling.util;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;

import java.awt.geom.Point2D;

/**
 * Created by yunyu on 10/2/15.
 */
public class FindMaximumForEnvelopeUsingBrent {
    double[] _branchLengthBound = {0.001, 6};
    double[] _inheritanceProbBound = {0.001, 0.999};
    boolean _printDetails = false;

    public void setBounds(double[] branchLengths, double[] inheritanceProb){
        _branchLengthBound = branchLengths;
        _inheritanceProbBound = inheritanceProb;
    }

    public Point2D getMaximum(Tuple<NetNode, NetNode> edge, Func<Double> likelihoodCalculator, Func2<Double,Boolean,Double> priorCalculator, double totalLnPriorMinus){
        double[] bounds = edge.Item1!=null?_branchLengthBound:_inheritanceProbBound;
        BrentOptimizer optimizer = new BrentOptimizer(0.000001, 0.000001);
        UnivariatePointValuePair result = optimizer.optimize(100, new FunctionToOptimize(edge, likelihoodCalculator, priorCalculator, totalLnPriorMinus), GoalType.MAXIMIZE, bounds[0], bounds[1]);
        if(_printDetails){
            System.out.println("Envelope: " + result.getPoint() + " resulting " + result.getValue());
            System.out.println(optimizer.getEvaluations());
        }
        return new Point2D.Double(result.getPoint(), result.getValue());
    }

    class FunctionToOptimize implements UnivariateFunction{
        Tuple<NetNode, NetNode> _edge;
        Func<Double> _likelihoodCalculator;
        Func2<Double,Boolean,Double> _priorCalculator;
        double _totalLnPriorMinus;

        public FunctionToOptimize(Tuple<NetNode, NetNode> edge, Func<Double> likelihoodCalculator, Func2<Double,Boolean,Double> priorCalculator, double totalLnPriorMinus) {
            _edge = edge;
            _likelihoodCalculator = likelihoodCalculator;
            _priorCalculator = priorCalculator;
            _totalLnPriorMinus = totalLnPriorMinus;

        }

        public double value(double sampleValue) {
            if (_edge.Item1 != null) {
                _edge.Item2.setParentDistance(_edge.Item1, sampleValue);
            } else {
                for (Object parentNode : _edge.Item2.getParents()) {
                    _edge.Item2.setParentProbability((NetNode) parentNode, sampleValue);
                    sampleValue = 1 - sampleValue;
                }
            }
            double lnLikelihood = _likelihoodCalculator.execute();
            double valueToMaximize = lnLikelihood + _totalLnPriorMinus + _priorCalculator.execute(sampleValue, _edge.Item1!=null);
            //System.out.println(sampleValue + " " + valueToMaximize);
            return valueToMaximize;
        }
    }

}
