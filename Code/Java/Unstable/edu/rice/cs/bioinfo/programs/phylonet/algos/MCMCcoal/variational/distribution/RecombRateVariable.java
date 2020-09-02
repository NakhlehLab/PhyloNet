package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution;


import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.RecombinationRate;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Variational distribution on recombination rate variable.
 *
 * Created by Xinhao Liu on 3/14/20.
 */
public class RecombRateVariable extends VariationalVariable {

    public RecombinationRate recombRate;

    public RecombRateVariable(double mean, double standardDeviation, RecombinationRate recombRate) {
        super(mean, standardDeviation);
        this.recombRate = recombRate;
    }

    @Override
    public void setVariableValue(double value) {
//        recombRate.setRecombRate(value * Utils.RECOMB_RATE_SCALE);
    }

    @Override
    public void meanGradientUpdate(double gradient) {
        setMean(getMean() + Utils.RECOMB_RATE_MEAN_LEARNING_RATE * gradient);
    }

    @Override
    public void standardDeviationGradientUpdate(double gradient) {
        setStandardDeviation(getStandardDeviation() + Utils.RECOMB_RATE_STDDEV_LEARNING_RATE * gradient);
    }
}
