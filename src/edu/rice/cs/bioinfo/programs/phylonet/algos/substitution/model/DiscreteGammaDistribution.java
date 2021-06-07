package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.special.Gamma;

import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 10/13/17
 * Time: 12:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class DiscreteGammaDistribution {
    private Random _random = null;
    private GammaDistribution _gamma = null;
    private double _shape = 1.0;
    private double _scale = 1.0;
    private int _num = 2;
    private double _freq[] = null;
    private double _rate[] = null;

    // use mean
    public DiscreteGammaDistribution(double shape, double scale, int num) {
        _shape = shape;
        _scale = scale;
        _num = num;
        _gamma = new GammaDistribution(shape, scale);
        _random = new Random();
        _rate = new double[_num];
        _freq = new double[_num];

        double cutPoints[] = new double[_num];
        cutPoints[0] = 0.0;
        for(int i = 1 ; i < _num ; i++) {
            cutPoints[i] = _gamma.inverseCumulativeProbability(1.0 * i / _num);
        }

        for(int i = 0 ; i < _num ; i++) {
            if(i == num - 1) {
                double I2 = Gamma.regularizedGammaP(_shape + 1, cutPoints[i] / _scale);
                _rate[i] = _shape * _scale * (1.0 - I2) * _num;
            } else {
                double I1 = Gamma.regularizedGammaP(_shape + 1, cutPoints[i + 1] / _scale);
                double I2 = Gamma.regularizedGammaP(_shape + 1, cutPoints[i] / _scale);
                _rate[i] = _shape * _scale * (I1 - I2) * _num;
            }
        }

        for(int i = 0 ; i < _num ; i++) {
            _freq[i] = 1.0 / _num;
        }

        for(int i = 1 ; i < _num ; i++) {
            _freq[i] += _freq[i - 1];
        }

        double sum = 0.0;
        for(int i = 0 ; i < _num ; i++) {
            sum += _rate[i] / _num;
        }

    }

    public void setSeed(Long seed) {
        _random = new Random(seed);
    }

    public double sample() {
        double r = _random.nextDouble();
        for(int i = 0 ; i < _num ; i++) {
            if(r < _freq[i])
                return _rate[i];
        }
        return 0.0;
    }

    public static void main(String[] Args)
    {
        DiscreteGammaDistribution dist = new DiscreteGammaDistribution(5, 0.2, 3);
        double sum = 0;
        for(int i = 0 ; i < 100000 ; i++) {
            sum += dist.sample();
        }
        System.out.println(sum);
    }
}
