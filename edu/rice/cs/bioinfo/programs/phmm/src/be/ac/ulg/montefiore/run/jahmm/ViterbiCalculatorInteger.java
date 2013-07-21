package be.ac.ulg.montefiore.run.jahmm;

import java.util.*;

import phylogeny.EvoTree;

import reader.*;

// kliu - pull in additional library support
import edu.rice.cs.bioinfo.library.programming.Tuple;


/* 
 * This class is the Viterbi Calculator that especially works for our integer observations.
 */

public class ViterbiCalculatorInteger {
	
	/*
	 * The psy and delta values, as described in Rabiner and Juand classical
	 * papers.
	 */
	private double[][] delta; 
	private int[][] psy;
	private int[] stateSequence;
	private double lnProbability;
	
	
	/**
	 * Computes the most likely state sequence matching an observation
	 * sequence given an HMM.
	 *
	 * @param hmm A Hidden Markov Model;
	 * @param oseq An ObservationInteger sequence;
	 * 
	 * @param treeStates All the trees(states) in this Hmm;
	 * @param seqType A map from species (of sequences) to integers;
	 * @param letterMap A map from alphabet letters to integers;
	 */
	public <O extends Observation> 
	ViterbiCalculatorInteger(List<ObservationInteger> oseq, 
				 ArrayList<Tuple<EvoTree,EvoTree>> treeStates, HashMap<String, Integer> seqType, HashMap<String, Integer> letterMap,
			Hmm<ObservationInteger> hmm)
	{
		if (oseq.isEmpty())
			throw new IllegalArgumentException("Invalid empty sequence");
		
		delta = new double[oseq.size()][hmm.nbStates()];
		psy = new int[oseq.size()][hmm.nbStates()];
		stateSequence = new int[oseq.size()];
		
		for (int i = 0; i < hmm.nbStates(); i++) {
			if (!checkEmission(hmm, i, new ObservationInteger(0))) {
				String tempStr = IntTranslator.intToLetter(seqType.size(), oseq.get(0).value, letterMap);
				Parser.inputEmissions(tempStr, 0, treeStates, seqType, hmm);
			}
			delta[0][i] = -Math.log(hmm.getPi(i)) - 
			Math.log(hmm.getOpdf(i).probability(oseq.get(0)));
			psy[0][i] = 0;
		}
		
		Iterator<ObservationInteger> oseqIterator = oseq.iterator();
		if (oseqIterator.hasNext())
			oseqIterator.next();
		
		int t = 1;
		
		while (oseqIterator.hasNext()) {
			ObservationInteger observation = oseqIterator.next();
			
			for (int i = 0; i < hmm.nbStates(); i++) {
				if (!checkEmission(hmm, i, observation)) {
					String tempStr = IntTranslator.intToLetter(seqType.size(), observation.value, letterMap);
					Parser.inputEmissions(tempStr, observation.value, treeStates, seqType, hmm);
				}
				computeStep(hmm, observation, t, i);
			}
			
			t++;
		}
		
		lnProbability = Double.MAX_VALUE;
		for (int i = 0; i < hmm.nbStates(); i++) {
			double thisProbability = delta[oseq.size()-1][i];
			
			if (lnProbability > thisProbability) {
				lnProbability = thisProbability;
				stateSequence[oseq.size() - 1] = i;
			}
		}
		lnProbability = -lnProbability;
		
		for (int t2 = oseq.size() - 2; t2 >= 0; t2--)
			stateSequence[t2] = psy[t2+1][stateSequence[t2+1]];
	}
	
	
	/*
	 * Computes delta and psy[t][j] (t > 0) 
	 */
	private <O extends Observation> void
	computeStep(Hmm<O> hmm, O o, int t, int j) 
	{
		double minDelta = Double.MAX_VALUE;
		int min_psy = 0;
		
		for (int i = 0; i < hmm.nbStates(); i++) {
			double thisDelta = delta[t-1][i] - Math.log(hmm.getAij(i, j));
			
			if (minDelta > thisDelta) {
				minDelta = thisDelta;
				min_psy = i;
			}
		}
		
		delta[t][j] = minDelta - Math.log(hmm.getOpdf(j).probability(o));
		psy[t][j] = min_psy;
	}
	
	
	/**
	 * Returns the neperian logarithm of the probability of the given
	 * observation sequence on the most likely state sequence of the given
	 * HMM.
	 *
	 * @return <code>ln(P[O,S|H])</code> where <code>O</code> is the given
	 *         observation sequence, <code>H</code> the given HMM and 
	 *         <code>S</code> the most likely state sequence of this observation
	 *         sequence given this HMM.
	 */
	public double lnProbability()
	{
		return lnProbability;
	}
	
	
	/**
	 * Returns a (clone of) the array containing the computed most likely
	 * state sequence.
	 *
	 * @return The state sequence; the i-th value of the array is the index
	 *         of the i-th state of the state sequence.
	 */
	public int[] stateSequence() 
	{
		return stateSequence.clone();
	}
	
	
	public static boolean checkEmission(Hmm<ObservationInteger> hmm, int stateNo, ObservationInteger obsNo) {
		OpdfInteger tempOpdf = (OpdfInteger) hmm.getOpdf(stateNo);
		if (tempOpdf.probability(obsNo) > 1) {
			// System.out.println("The obsNo is : " + obsNo);
			return false;
		}
		else
			return true;
	}

}
