package edu.rice.cs.bioinfo.programs.phylonet.algos.gibbssampling.util;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import org.freehep.math.minuit.FCNBase;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMinimize;
import org.freehep.math.minuit.MnUserParameters;

import java.awt.geom.Point2D;

/**
 * Created by yunyu on 10/2/15.
 */
public class FindMaximumForEnvelopeUsingMinuit {
    double[] _branchLengthBound = {0.001, 6};
    double[] _inheritanceProbBound = {0.000001, 0.999999};
    boolean _printDetails = false;

    public Point2D getMaximum(Tuple<NetNode, NetNode> edge, Func<Double> likelihoodCalculator, Func2<Double,Boolean,Double> priorCalculator, double totalLnPriorMinus){
        MnUserParameters upar = new MnUserParameters();
        if(edge.Item1!=null) {
            upar.add("x", edge.Item2.getParentDistance(edge.Item1), 1, _branchLengthBound[0], _branchLengthBound[1]);
            //upar.setLowerLimit(0, 0);
        }else{
            upar.add("x", edge.Item2.getParentProbability((NetNode) edge.Item2.getParents().iterator().next()), 1, _inheritanceProbBound[0], _inheritanceProbBound[1]);
        }
        FindMaximum fcn = new FindMaximum(edge, likelihoodCalculator, priorCalculator, totalLnPriorMinus);
        MnMinimize migrad = new MnMinimize(fcn, upar);
        //migrad.setPrecision(0.000001);
        FunctionMinimum min = migrad.minimize();
        double maximumValue = min.fval() * (-1);
        if(_printDetails){
            System.out.println("Envelope: " + min.userParameters().value(0) + " resulting " + maximumValue);
            //System.out.println(migrad.numOfCalls());
            //System.exit(-1);
        }
        return new Point2D.Double(min.userParameters().value(0), maximumValue);
    }

    public void setBounds(double[] branchLengths, double[] inheritanceProb){
        _branchLengthBound = branchLengths;
        _inheritanceProbBound = inheritanceProb;
    }

    class FindMaximum implements FCNBase {
        Tuple<NetNode, NetNode> _edge;
        Func<Double> _likelihoodCalculator;
        Func2<Double,Boolean,Double> _priorCalculator;
        double _totalLnPriorMinus;

        public FindMaximum(Tuple<NetNode, NetNode> edge, Func<Double> likelihoodCalculator, Func2<Double,Boolean,Double> priorCalculator, double totalLnPriorMinus) {
            _edge = edge;
            _likelihoodCalculator = likelihoodCalculator;
            _priorCalculator = priorCalculator;
            _totalLnPriorMinus = totalLnPriorMinus;

        }

        public double valueOf(double[] par) {
            double sampleValue = par[0];
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
            //System.out.println(sampleValue + ": " + valueToMaximize);
            return (-1) * valueToMaximize;
        }




    }
}
