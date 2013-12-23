/**
 * Simple container class for transition probability parameters
 * other than basic coalescent model parameters (i.e., parental
 * tree branch lengths).
 *
 * See writeup for details.
 */

package be.ac.ulg.montefiore.run.jahmm.phmm;

public class TransitionProbabilityParameters {
    // RECOMBINATION_FREQUENCY
    public enum ParameterChoice { HYBRIDIZATION_FREQUENCY }
    //protected double recombinationFrequency;
    protected double hybridizationFrequency;
    
    public TransitionProbabilityParameters (double inHybridizationFrequency) {
	setHybridizationFrequency(inHybridizationFrequency);
    }
    
    public double get (ParameterChoice parameterChoice) {
	double value;
	switch (parameterChoice) {
	// case RECOMBINATION_FREQUENCY:
	//     value = recombinationFrequency;
	//     break;
	case HYBRIDIZATION_FREQUENCY:
	default:
	    value = hybridizationFrequency;
	    break;
	}
	
	return (value);
    }

    public void set (ParameterChoice parameterChoice,
		     double inValue) {
	switch (parameterChoice) {
	// case RECOMBINATION_FREQUENCY:
	//     recombinationFrequency = inValue;
	//     break;
	case HYBRIDIZATION_FREQUENCY:
	default:
	    hybridizationFrequency = inValue;
	    break;
	}
    }

    public double getHybridizationFrequency () {
	return (get(ParameterChoice.HYBRIDIZATION_FREQUENCY));
    }

    public void setHybridizationFrequency (double inHybridizationFrequency) {
	    set(ParameterChoice.HYBRIDIZATION_FREQUENCY, inHybridizationFrequency);
    }
}

    // double inRecombinationFrequency,
    //setRecombinationFrequency(inRecombinationFrequency);

    // public double getRecombinationFrequency () {
    // 	return (get(ParameterChoice.RECOMBINATION_FREQUENCY));
    // }

    // public void setRecombinationFrequency (double inRecombinationFrequency) {
    // 	set(ParameterChoice.RECOMBINATION_FREQUENCY, inRecombinationFrequency);
    // }
    
