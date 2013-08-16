package optimize;

/**
 * Directly handle model manipulation in here.
 */

import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;
import runHmm.runHmm;

public class SwitchingFrequencyParameter extends Parameter {
    public static final double DEFAULT_MINIMUM_PROBABILITY = MultivariateOptimizer.DEFAULT_MINIMUM_PROBABILITY;
    public static final double DEFAULT_INITIAL_PROBABILITY = MultivariateOptimizer.DEFAULT_INITIAL_PROBABILITY;
    // need to calculate it on the basis of the current HMM parameter values
    public static final double DEFAULT_MAXIMUM_PROBABILITY = MultivariateOptimizer.DEFAULT_MAXIMUM_PROBABILITY;

    protected runHmm runHmmObject; // for pushing updated transition probabilities
    protected TransitionProbabilityParameters transitionProbabilityParameters;
    protected TransitionProbabilityParameters.ParameterChoice parameterChoice;

    // need to access rest of HMM state to figure out maximum possible frequency parameter value
    //protected runHmm runHmmObject;

    /**
     * User must explicitly request whether or not an update is requested 
     * during construction.
     *
     * Don't really need access to HMM state anymore.
     */
    public SwitchingFrequencyParameter (String inName, 
					double inValue, 
					runHmm inRunHmmObject,
					TransitionProbabilityParameters inTransitionProbabilityParameters, 
					TransitionProbabilityParameters.ParameterChoice inParameterChoice, 
					boolean checkValueMinimumConstraintFlag,
					boolean checkValueMaximumConstraintFlag,
					boolean updateModelStateFlag) {
	// order forced by Java language constraints
	// need to delay the setValue() check until after runHmmObject reference ready
	super(inName, inValue, checkValueMinimumConstraintFlag, checkValueMinimumConstraintFlag, false);

	this.runHmmObject = inRunHmmObject;
	this.transitionProbabilityParameters = inTransitionProbabilityParameters;
	this.parameterChoice = inParameterChoice;

	if (updateModelStateFlag) {
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
	// need to refer to current HMM parameter value settings to figure 
	// out maximum possible frequency parameter setting
	return (DEFAULT_MAXIMUM_PROBABILITY);
    }

    public void updateModelState () {
	transitionProbabilityParameters.set(parameterChoice, getValue());

	// propagate hybridization/etc. frequency parameters on to
	// HMM transition probability matrix
	runHmmObject.updateTransitionProbabilities();
    }

    public TransitionProbabilityParameters.ParameterChoice getParameterChoice () {
	return (parameterChoice);
    }


}

// runHmm inRunHmmObject, 

	//this.runHmmObject = inRunHmmObject;

	// retry the setValue now
	// since runHmmObject reference finally ready
	//setValue(inValue, updateFlag);

	//return (runHmmObject.calculateMaximumFrequencyParameter(getParameterChoice()));

    /**
     * WARNING: doesn't update attached TransitionProbabilityParameters object by default!
     * Issue with constructor and super().
     */
    // public void setValue (double inValue) {
    // 	this.setValue(inValue, false);
    // }

    // public void setValue (double inValue, boolean updateFlag) {
    // 	super.setValue(inValue);
    // 	if (updateFlag) {
    // 	    updateModelState();
    // 	}
    // }
