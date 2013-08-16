package optimize;

import substitutionModel.GTRSubstitutionModel;
import java.util.Vector;

public class GTRBaseFrequencyParameter extends Parameter {
    protected GTRSubstitutionModel gtrSubstitutionModel;
    // choice of which parameter to optimize
    // use relative weights for all
    // by convention - the first is always set to 1.0
    protected int index;
    // need to keep track of all base frequency parameters
    // weighting of one is relative to weighting of others
    //
    // the first is always set to 1.0
    protected Vector<GTRBaseFrequencyParameter> gbfps;
    protected CalculationCache calculationCache;

    public GTRBaseFrequencyParameter (String inName, 
				      double inValue,
				      GTRSubstitutionModel inGTRSubstitutionModel,
				      int inIndex,
				      Vector<GTRBaseFrequencyParameter> inGbfps,
				      CalculationCache inCalculationCache,
				      boolean checkValueMinimumConstraintFlag,
				      boolean checkValueMaximumConstraintFlag,
				      boolean updateModelStateFlag) {
	super(inName, inValue, checkValueMinimumConstraintFlag, checkValueMaximumConstraintFlag, false);

	this.gtrSubstitutionModel = inGTRSubstitutionModel;
	// paranoid
	if ((inIndex < 0) || (inIndex >= gtrSubstitutionModel.getAlphabet().length())) {
	    throw (new RuntimeException("ERROR: index in GTRRateParameter(...) is out of bounds. " + inIndex));
	}
	this.index = inIndex;
	this.gbfps = inGbfps;
	this.calculationCache = inCalculationCache;

	if (updateModelStateFlag) {
	    updateModelState();
	}
    }

    public void updateModelState () {
	// all base frequencies must change as well
	double[] newFreqParameters = calculateValueOfAllGTRBaseFrequencyParameters();
	double[] originalFreqParameters = gtrSubstitutionModel.getStationaryProbabilities();
	// paranoid
	if (newFreqParameters.length != originalFreqParameters.length) {
	    throw (new RuntimeException("ERROR: newFreqParameters.length is not equal to originalFreqParameters.length in updateModelState()."));
	}
	// set by reference
	for (int i = 0; i < originalFreqParameters.length; i++) {
	    originalFreqParameters[i] = newFreqParameters[i];
	}
	gtrSubstitutionModel.updateRateMatrix();

	// totally clear out associated caches
	calculationCache.cacheSubstitutionProbabilityMatrix.clear();
	calculationCache.cacheSubstitutionProbability.clear();
    }

    /**
     * Calculates base frequencies based on current set of parameter weights.
     */
    protected double[] calculateValueOfAllGTRBaseFrequencyParameters () {
	double totalWeight = 0.0;
	for (GTRBaseFrequencyParameter gbfp : gbfps) {
	    totalWeight += gbfp.getValue();
	}
	// solve for x
	double x = 1.0 / totalWeight;
	double[] result = new double[gbfps.size()];
	for (int i = 0; i < gbfps.size(); i++) {
	    result[i] = x * gbfps.get(i).getValue();
	}
	return (result);
    }

    /**
     * Since all we need are weights for these parameters, just overload
     * rate default settings.
     */
    public double getMinimumValue () {
	return (MultivariateOptimizer.DEFAULT_MINIMUM_RATE);
    }

    public double getDefaultInitialValue () {
	return (MultivariateOptimizer.DEFAULT_INITIAL_RATE);
    }

    public double getMaximumValue () {
	return (MultivariateOptimizer.DEFAULT_MAXIMUM_RATE);
    }


}
