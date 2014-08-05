package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util;

import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model.JahmmNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;

import java.util.List;

public class StateProbabilityCalculator extends BaumWelchScaledLearner {

    /**
     * This calculates a two dimensional array storing the probability of being in a state at time t.
     * @param obsSeq The observations.
     * @param hmm The hidden markov model.
     * @return A two dimensional array, where arr[time][state] is the probability of being in state at time.
     */
    public double[][] getProbabilityOfStateAtTime(List<JahmmNucleotideObservation> obsSeq, Hmm<JahmmNucleotideObservation> hmm)
    {
        ForwardBackwardCalculator fbc =
                generateForwardBackwardCalculator(obsSeq, hmm);

        double xi[][][] = estimateXi(obsSeq, fbc, hmm);
        return  estimateGamma(xi, fbc);
    }
}
