package gridSearch;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

import phylogeny.EvoTree;
import runHmm.AllInformation;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;

/**
 * This knob is used for the recombination frequency parameter from
 * the TransitionProbabilityParameters class.
 * @author k3kathy
 *
 */
public class RecombinationFreqNob extends Nob {

	private Hmm thisHmm;
	
	private TransitionProbabilityParameters probsParam;
	private ArrayList<HiddenState> trees_states;
	private Map<EvoTree,Set<HiddenState>> parentalTreeClasses;
	
	public RecombinationFreqNob(int gIn, double minIn, double maxIn, 
			Hmm hmmIn, TransitionProbabilityParameters probsParamIn, ArrayList<HiddenState> treeStatesIn, 
			Map<EvoTree, Set<HiddenState>> parentalTreeIn) {
		super(gIn, minIn, maxIn);
		this.thisHmm = hmmIn;
		this.probsParam = probsParamIn;
		this.trees_states = treeStatesIn;
		this.parentalTreeClasses = parentalTreeIn;
	}
	
	public void set_param(double value) {
		probsParam.setRecombinationFrequency(value);
		double[][] newTransition = AllInformation.calculateAij(trees_states, value, probsParam.getHybridizationFrequency(),
				parentalTreeClasses);
		thisHmm.setTransitionMatrix(newTransition);
	}
	

	public double get_param() {
		return probsParam.getRecombinationFrequency();
	}
	
}
