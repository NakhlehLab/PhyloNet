package optimize;

import substitutionModel.GTRSubstitutionModel;

public class GTRRateParameter extends Parameter {
    public static final double DEFAULT_MINIMUM_RATE = MultivariateOptimizer.DEFAULT_MINIMUM_RATE;
    public static final double DEFAULT_INITIAL_RATE = MultivariateOptimizer.DEFAULT_INITIAL_RATE;
    public static final double DEFAULT_MAXIMUM_RATE = MultivariateOptimizer.DEFAULT_MAXIMUM_RATE;

    protected GTRSubstitutionModel gtrSubstitutionModel;
    // choice of which parameter to optimize
    protected int index;

    public GTRRateParameter (String inName, 
			     double inValue,
			     GTRSubstitutionModel inGTRSubstitutionModel,
			     int inIndex,
			     boolean checkValueMinimumConstraintFlag,
			     boolean checkValueMaximumConstraintFlag,
			     boolean updateModelStateFlag) {
	super(inName, inValue, checkValueMinimumConstraintFlag, checkValueMaximumConstraintFlag, false);

	this.gtrSubstitutionModel = inGTRSubstitutionModel;
	// paranoid
	if ((inIndex < 0) || (inIndex >= gtrSubstitutionModel.getRateParameterCount())) {
	    throw (new RuntimeException("ERROR: index in GTRRateParameter(...) is out of bounds. " + inIndex));
	}
	this.index = inIndex;
	
	if (updateModelStateFlag) {
	    updateModelState();
	}
    }

    public double getMinimumValue () {
	return (DEFAULT_MINIMUM_RATE);
    }

    public double getDefaultInitialValue () {
	return (DEFAULT_INITIAL_RATE);
    }

    public double getMaximumValue () {
	return (DEFAULT_MAXIMUM_RATE);
    }

    public void updateModelState () {
	double[] originalRateParameters = gtrSubstitutionModel.getOriginalRateParameters();
	// set by reference
	originalRateParameters[index] = this.getValue();
	gtrSubstitutionModel.updateRateMatrix();
    }
}
