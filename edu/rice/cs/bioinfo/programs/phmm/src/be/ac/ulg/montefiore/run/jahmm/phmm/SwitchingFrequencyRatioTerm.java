/**
 * Simple container class for ratio term (positive real number)
 * corresponding to a switching frequency parameter used to calculate
 * hidden state transition probabilities.
 *
 * No storage of which trees (parental tree or gene genealogy)
 * this corresponds to. Maintained externally using a BidirectionalMultimap.
 *
 * See writeup for details.
 */

package be.ac.ulg.montefiore.run.jahmm.phmm;

import optimize.CalculationCache;

public class SwitchingFrequencyRatioTerm {
    protected String name;
    protected double value;
    // perform cache invalidation in here
    protected CalculationCache calculationCache;
    // set to true to invalidate cacheParentalTreeSwitchingFrequencyMap
    // otherwise invalidates cacheGeneGenealogySwitchingFrequencyMap
    protected boolean invalidateParentalTreeSwitchingFrequencyMapFlag;
    // bleh, also need to store the maximum number of alternatives
    // (parental trees, gene genealogies, or whatever)
    // to help cap the switching frequencies
    protected int numAlternatives;
    
    public SwitchingFrequencyRatioTerm (String inName, double inValue, CalculationCache inCalculationCache, boolean inInvalidateParentalTreeSwitchingFrequencyMapFlag, int inNumAlternatives) {
	this.calculationCache = inCalculationCache;
	this.invalidateParentalTreeSwitchingFrequencyMapFlag = inInvalidateParentalTreeSwitchingFrequencyMapFlag;
	this.numAlternatives = inNumAlternatives;

	setName(inName);
	setValue(inValue);
    }
    
    public String getName () {
	return (name);
    }
    
    public double getValue () {
	return (value);
    }

    public void setName (String inName) {
	name = inName;
    }

    public void setValue (double inValue) {
	if (inValue > 0.0) {
	    value = inValue;
	    // wipe cache
	    if (invalidateParentalTreeSwitchingFrequencyMapFlag) {
		calculationCache.cacheParentalTreeSwitchingFrequencyMap = null;
	    }
	    else {
		calculationCache.cacheGeneGenealogySwitchingFrequencyMap = null;
	    }
	}
	else {
	    throw (new RuntimeException("ERROR: SwitchingFrequencyRatioTerm.setValue(double) called with non-positive number."));
	}
    }

    public int getNumAlternatives () {
	return (numAlternatives);
    }

    /**
     * For map support.
     */
    public boolean equals (Object obj) {
	if (!(obj instanceof SwitchingFrequencyRatioTerm)) {
	    return (false);
	}

	SwitchingFrequencyRatioTerm so = (SwitchingFrequencyRatioTerm) obj;
	return (this.getName().equals(so.getName()));
    }

    public int hashCode() {
	return (this.getName().hashCode());
    }

    public String toString () {
	return (getName() + ": " + getValue());
    }
}
