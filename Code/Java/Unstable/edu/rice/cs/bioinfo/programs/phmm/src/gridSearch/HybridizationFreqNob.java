package gridSearch;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

import phylogeny.EvoTree;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;

/**
 * This knob is used for the hybridization frequency parameter from
 * the TransitionProbabilityParameters class.
 * @author k3kathy
 *
 */
public class HybridizationFreqNob extends Nob {

	private Hmm thisHmm;
	
	private TransitionProbabilityParameters probsParam;
	private ArrayList<HiddenState> trees_states;
	private Map<EvoTree,Set<HiddenState>> parentalTreeClasses;
	
	public HybridizationFreqNob(int gIn, double minIn, double maxIn, TransitionProbabilityParameters probsParamIn) {
		super(gIn, minIn, maxIn);
		
		this.probsParam = probsParamIn;
	}
	
	public void set_param(double value) {
		probsParam.setHybridizationFrequency(value);
		double[][] newTransition = GridSearchAlgorithm.calculateAij(trees_states, probsParam.getRecombinationFrequency(), value,
				parentalTreeClasses);
		thisHmm.setTransitionMatrix(newTransition);
	}
	

	public double get_param() {
		return probsParam.getHybridizationFrequency();
	}
	
}
