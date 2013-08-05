package gridSearch;

import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;

/**
 * This knob is used for the recombination frequency parameter from
 * the TransitionProbabilityParameters class.
 * @author k3kathy
 *
 */
public class RecombinationFreqNob extends Nob {

	private TransitionProbabilityParameters probsParam;
	
	public RecombinationFreqNob(int gIn, double minIn, double maxIn, TransitionProbabilityParameters probsParamIn) {
		super(gIn, minIn, maxIn);
		
		this.probsParam = probsParamIn;
	}
	
	public void set_param(double value) {
		probsParam.setRecombinationFrequency(value);
	}
	

	public double get_param() {
		return probsParam.getRecombinationFrequency();
	}
	
}
