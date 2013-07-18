package be.ac.ulg.montefiore.run.jahmm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import phylogeny.EvoTree;

public class MyHMM {
	
	public static Hmm<ObservationInteger> buildMyHmm(double[] pi, double[][] a, int nbObs) {
		
		if (pi.length != a.length)
			throw new IllegalArgumentException("The dimension of the initial probability array " +
					"does not match the dimension of the transition probability matrix.");
		
		int nbStates = pi.length;
		
		ArrayList<Opdf<ObservationInteger>> opdfs = new ArrayList<Opdf<ObservationInteger>>(nbStates);
		
		for (int i = 0; i < nbStates; i++) {
			double[] probabilities = new double[nbObs];
			for (int j = 0; j < nbObs; j++) {
				probabilities[j] = 10;
			}
			Opdf<ObservationInteger> temp = new OpdfInteger(probabilities);
			opdfs.add(temp);
		}
		
		Hmm<ObservationInteger> newHmm = new Hmm<ObservationInteger>(pi, a, opdfs);
		
		return newHmm;
		
	}
	
	public static void setEmission(Hmm<ObservationInteger> curHmm, int state, int obsNo, double prob) {
		OpdfInteger curState = (OpdfInteger)curHmm.getOpdf(state);
		curState.setProb(obsNo, prob);
		
	}
	
	
	/**
	 * Returns an array containing the most likely state sequence matching an
	 * observation sequence given this HMM.  This sequence <code>I</code>
	 * maximizes the probability of <code>P[I|O,Model]</code> where
	 * <code>O</code> is the observation sequence and <code>Model</code> this
	 * HMM model.
	 *
	 * @param oseq A non-empty observation sequence.
	 * @return An array containing the most likely sequence of state numbers.
	 *         This array can be modified.
	 */
	public static int[] mostLikelyStateSequenceForIntObs(List<ObservationInteger> oseq, 
			ArrayList<EvoTree> treeStates, HashMap<String, Integer> seqType, HashMap<String, Integer> letterMap,
			Hmm<ObservationInteger> hmm)
	{
		return (new ViterbiCalculatorInteger(oseq, treeStates, seqType, letterMap, hmm)).stateSequence();
	}

}
