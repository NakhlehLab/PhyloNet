/**
 * Simple container class for ratio term (positive real number)
 * corresponding to a switching frequency parameter used to calculate
 * hidden state transition probabilities.
 *
 * No storage of which trees (parental tree or gene genealogy)
 * this corresponds to. Maintained externally using a BidirectionalMap.
 *
 * See writeup for details.
 */

package be.ac.ulg.montefiore.run.jahmm.phmm;

public class SwitchingFrequencyRatioTerm {
    protected double term;
    
    public SwitchingFrequencyRatioTerm (double inTerm) {
	set(inTerm);
    }
    
    public double get () {
	return (term);
    }

    public void set (double inTerm) {
	if (inTerm > 0.0) {
	    term = inTerm;
	}
	else {
	    throw (new RuntimeException("ERROR: SwitchingFrequencyRatioTerm.term(double) called with non-positive number."));
	}
    }
}
