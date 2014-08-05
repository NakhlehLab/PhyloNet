package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;

import java.util.HashMap;
import java.util.Map;


public abstract class NucleotideProbabilityAlgorithm {
    public abstract double getProbability(NucleotideObservation dna);

    public Map<NucleotideObservation,Double> memo = new HashMap<NucleotideObservation,Double>();

    public double getProbabilityCached(NucleotideObservation dna)
    {
        Double result = memo.get(dna);
        if (result == null)
        {
            result = getProbability(dna);
            memo.put(dna,result);
            return result;
        }
        return result;
    }

    public void seedProbability(NucleotideObservation dna, double probability)
    {
        memo.put(dna,probability);
    }
}
