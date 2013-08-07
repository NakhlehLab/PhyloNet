package optimize;

/**
 * Directly handle model manipulation in here.
 */

import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;

public class FrequencyParameter extends Parameter {
    public static final double DEFAULT_MINIMUM_PROBABILITY = MultivariateOptimizer.DEFAULT_MINIMUM_PROBABILITY;
    public static final double DEFAULT_INITIAL_PROBABILITY = MultivariateOptimizer.DEFAULT_INITIAL_PROBABILITY;
    public static final double DEFAULT_MAXIMUM_PROBABILITY = MultivariateOptimizer.DEFAULT_MAXIMUM_PROBABILITY;

    protected TransitionProbabilityParameters.ParameterChoice parameterChoice;

    protected TransitionProbabilityParameters transitionProbabilityParameters;

    /**
     * User must explicitly request whether or not an update is requested 
     * during construction.
     */
    public FrequencyParameter (String inName, double inValue, TransitionProbabilityParameters inTransitionProbabilityParameters, TransitionProbabilityParameters.ParameterChoice inParameterChoice, boolean updateFlag) {
	// order forced by Java language constraints
    super(inName, inValue);

	this.parameterChoice = inParameterChoice;
	this.transitionProbabilityParameters = inTransitionProbabilityParameters;

	if (updateFlag) {
	    updateModelState();
	}
    }

    public double getMinimumValue () {
    return (DEFAULT_MINIMUM_PROBABILITY);
    }

    public double getDefaultInitialValue () {
	return (DEFAULT_INITIAL_PROBABILITY);
    }

    public double getMaximumValue () {
    return (DEFAULT_MAXIMUM_PROBABILITY);
    }

    /**
     * WARNING: doesn't update attached TransitionProbabilityParameters object by default!
     * Issue with constructor and super().
     */
    public void setValue (double inValue) {
	this.setValue(inValue, false);
    }

    public void setValue (double inValue, boolean updateFlag) {
	super.setValue(inValue);
	if (updateFlag) {
	    updateModelState();
	}
    }

    public void updateModelState () {
	transitionProbabilityParameters.set(parameterChoice, getValue());
    }

    public TransitionProbabilityParameters.ParameterChoice getParameterChoice () {
	return (parameterChoice);
    }


}
