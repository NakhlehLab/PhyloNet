package optimize;

/**
 * Directly handle model manipulation in here.
 * Lightweight. Fiddles with SwitchingFrequencyRatio as a simple positive rate (member of a ratio).
 */

import be.ac.ulg.montefiore.run.jahmm.phmm.SwitchingFrequencyRatioTerm;
import runHmm.runHmm;

public class SwitchingFrequencyRatioTermParameter extends Parameter {
    // allow this to get basically to zero
    // divide the max by 10^10
    //public static final double DEFAULT_MINIMUM_RATE = MultivariateOptimizer.ABSOLUTE_ACCURACY;
    // careful - max depends on SwitchingFrequencyRatioTerm.getNumAlternatives()
    // eh - just go an order of magnitude below max
    //public static final double DEFAULT_INITIAL_RATE = MultivariateOptimizer.DEFAULT_INITIAL_RATE;
    // Like in original algorithm,
    // switching frequency approaching half or more can cause a degenerate corner case.
    // Lots of switching artifacts.
    //
    // Cap this so that maximum switching frequency is 0.25 .
    //
    // For two parental tree case, term is maximum  1/3 .
    //
    // Without any complicated structure to the ratio, the maximum is 1/(3(k-1)) where k is 
    // the number of alternative "states" (number of parental trees, or gene genealogies, or whatever).
    //
    // Conservatively set any non-self-transition probability to be at most 1/(3(k-1)),
    // ensuring that the self-transition probability is at least 3/4.
    //public static final double DEFAULT_MAXIMUM_RATE = MultivariateOptimizer.DEFAULT_MAXIMUM_RATE;

    // corresponds to maximum total non-self-transition switching frequency of 1/3
    //
    // actual constraint is applied in runHmm.calculateSwitchingFrequencies(...)
    //public static final double DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_TERMS_TOTAL_WEIGHT = 1.0 / 2.0;

    // rest calculated off of maximum
    public static final double DEFAULT_MINIMUM_RATIO_TERM = 1e-6;
    public static final double DEFAULT_INITIAL_RATIO_TERM = 1e-3;
    // any more than this will certainly cause normalization+rescaling in runHmm.calculateSwitchingFrequencies(...)
    //public static final double DEFAULT_MAXIMUM_RATIO_TERM = DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_TERMS_TOTAL_WEIGHT;
    public static final double DEFAULT_MAXIMUM_RATIO_TERM = MultivariateOptimizer.DEFAULT_MAXIMUM_RATE;

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

    // really, can go to zero
    public double getMinimumValue () {
	return (DEFAULT_MINIMUM_RATIO_TERM);
    }

    public double getDefaultInitialValue () {
	return (DEFAULT_INITIAL_RATIO_TERM);
    }

    public double getMaximumValue () {
	return (DEFAULT_MAXIMUM_RATIO_TERM);
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
