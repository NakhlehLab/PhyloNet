/**
 * Simple container class for transition probability parameters
 * other than basic coalescent model parameters (i.e., parental
 * tree branch lengths).
 *
 * See writeup for details.
 */

package be.ac.ulg.montefiore.run.jahmm.phmm;

public class TransitionProbabilityParameters {
    protected double recombinationFrequency;
    protected double hybridizationFrequency;

    public TransitionProbabilityParameters (double inRecombinationFrequency,
						  double inHybridizationFrequency) {
	setRecombinationFrequency(inRecombinationFrequency);
	setHybridizationFrequency(inHybridizationFrequency);
    }

    public double getRecombinationFrequency () {
	return (recombinationFrequency);
    }

    public void setRecombinationFrequency (double inRecombinationFrequency) {
	this.recombinationFrequency = inRecombinationFrequency;
    }
    
    public double getHybridizationFrequency () {
	return (hybridizationFrequency);
    }

    public void setHybridizationFrequency (double inHybridizationFrequency) {
	this.hybridizationFrequency = inHybridizationFrequency;
    }
}
