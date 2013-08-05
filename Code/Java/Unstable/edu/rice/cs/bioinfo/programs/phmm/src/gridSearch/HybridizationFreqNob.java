package gridSearch;

import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;


public class HybridizationFreqNob extends Nob {

	private TransitionProbabilityParameters probsParam;
	
	public HybridizationFreqNob(int gIn, double minIn, double maxIn, TransitionProbabilityParameters probsParamIn) {
		super(gIn, minIn, maxIn);
		
		this.probsParam = probsParamIn;
	}
	
	public void set_param(double value) {
		probsParam.setHybridizationFrequency(value);
	}
	

	public double get_param() {
		return probsParam.getHybridizationFrequency();
	}
	
}
