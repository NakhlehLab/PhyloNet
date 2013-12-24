package optimize;

/**
 * Directly handle model manipulation in here.
 * Lightweight. Fiddles with SwitchingFrequencyRatio as a simple positive rate (member of a ratio).
 */

import be.ac.ulg.montefiore.run.jahmm.phmm.SwitchingFrequencyRatioTerm;
import runHmm.runHmm;

public class SwitchingFrequencyRatioTermParameter extends Parameter {
    public static final double DEFAULT_MINIMUM_RATE = MultivariateOptimizer.DEFAULT_MINIMUM_RATE;
    public static final double DEFAULT_INITIAL_RATE = MultivariateOptimizer.DEFAULT_INITIAL_RATE;
    public static final double DEFAULT_MAXIMUM_RATE = MultivariateOptimizer.DEFAULT_MAXIMUM_RATE;

    protected runHmm runHmmObject; // for pushing updated transition probabilities
    protected SwitchingFrequencyRatioTerm sfrt;

    /**
     * User must explicitly request whether or not an update is requested 
     * during construction.
     */
    public SwitchingFrequencyRatioTermParameter (String inName, 
					double inValue, 
					runHmm inRunHmmObject,
					SwitchingFrequencyRatioTerm inSfrt,
					boolean checkValueMinimumConstraintFlag,
					boolean checkValueMaximumConstraintFlag,
					boolean updateModelStateFlag) {
	// order forced by Java language constraints
	// need to delay the setValue() check until after runHmmObject reference ready
	super(inName, inValue, checkValueMinimumConstraintFlag, checkValueMinimumConstraintFlag, false);

	this.runHmmObject = inRunHmmObject;
	this.sfrt = inSfrt;

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
	// invalidate caches
	// caches re-populated in updateTransitionProbabilities call below
	//
	// push to SwitchingFrequencyRatioTerm.setValue()
	sfrt.setValue(getValue());

	// propagate hybridization/etc. frequency parameters on to
	// HMM transition probability matrix
	//
	// Throw in caching into runHmm?
	// If parental branch lengths change, no need to 
	// re-calculate switching frequencies from switching-frequency-ratio-terms.
	runHmmObject.updateTransitionProbabilities();
    }

}
