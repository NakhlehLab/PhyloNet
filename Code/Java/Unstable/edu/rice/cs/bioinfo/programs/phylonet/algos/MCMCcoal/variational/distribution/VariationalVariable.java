package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution;


import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Variational distribution for one model parameter.
 * We use univariate normal distribution as variational distribution.
 *
 * Created by Xinhao Liu on 3/14/20.
 */
public abstract class VariationalVariable {

    protected NormalDistribution variationalDistribution;
    protected double mean;
    protected double standardDeviation;

    public VariationalVariable(double mean, double standardDeviation) {
        this.mean = mean;
        this.standardDeviation = standardDeviation;
        variationalDistribution = new NormalDistribution(this.mean, this.standardDeviation);
    }

    public double logDensity(double x) {
        return variationalDistribution.logDensity(x);
    }

    public double scoreMean(double x) {
        return (x - mean) / (standardDeviation * standardDeviation);
    }

    public double scoreStandardDeviation(double x) {
        return (-1.0 / standardDeviation) + (((x - mean) * (x - mean)) / (standardDeviation * standardDeviation * standardDeviation));
    }

    public double sample() {
        return variationalDistribution.sample();
    }

    public double getMean() {
        return mean;
    }

    public double getStandardDeviation() {
        return standardDeviation;
    }

    public void setMean(double newMean) {
        mean = newMean;
        variationalDistribution = new NormalDistribution(mean, standardDeviation);
    }

    public void setStandardDeviation(double newStandardDeviation) {
        standardDeviation = newStandardDeviation;
        variationalDistribution = new NormalDistribution(mean, standardDeviation);
    }

    public abstract void setVariableValue(double value);

    public abstract void meanGradientUpdate(double gradient);

    public abstract void standardDeviationGradientUpdate(double gradient);

}
