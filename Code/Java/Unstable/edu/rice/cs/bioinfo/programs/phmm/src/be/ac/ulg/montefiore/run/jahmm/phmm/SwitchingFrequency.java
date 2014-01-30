/**
 * Simple container class for switching frequency.
 *
 * See writeup for details.
 */

package be.ac.ulg.montefiore.run.jahmm.phmm;

import optimize.CalculationCache;

public class SwitchingFrequency extends SwitchingFrequencyRatioTerm {
    public static final String GAMMA = "gamma";

    public SwitchingFrequency (String inName, double inValue, CalculationCache inCalculationCache, int inNumAlternatives) {
	super(inName, inValue, inCalculationCache, inNumAlternatives);
    }

    public void setValue (double inValue) {
	if ((inValue >= 0.0) & (inValue <= 1.0)) {
	    value = inValue;
	    // wipe cache
	    //calculationCache.cacheSwitchingFrequencyMap = null;

	    // if (invalidateParentalTreeSwitchingFrequencyMapFlag) {
	    // 	calculationCache.cacheParentalTreeSwitchingFrequencyMap = null;
	    // }
	    // else {
	    // 	calculationCache.cacheGeneGenealogySwitchingFrequencyMap = null;
	    // }
	}
	else {
	    throw (new RuntimeException("ERROR: SwitchingFrequency.setValue(double) called with argument outside of [0,1] range."));
	}
    }

}

