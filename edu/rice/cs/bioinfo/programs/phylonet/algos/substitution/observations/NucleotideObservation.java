package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

public abstract class NucleotideObservation
{
    /**
     * Creates a list of NucletideObservations from a mapping of allele to the entire dna.
     * @param alleleToDNA A mapping of allele names to their associated dna.
     * @return A list of NucleotideObservations.
     */
    public static List<NucleotideObservation> getObservations(Map<String,String> alleleToDNA)
    {
        // I need the first dna object in order to determine the length.
        String firstDna = alleleToDNA.values().iterator().next();

        return getObservations(alleleToDNA,0,firstDna.length());
    }

    public static List<NucleotideObservation> getObservations(Map<String,String> alleleToDNA,int startIndex,int endIndex)
    {
        List<NucleotideObservation> result = new ArrayList<NucleotideObservation>();
        for (int i = startIndex; i<endIndex; i++)
            result.add(new ManyNucleotideObservation(alleleToDNA,i));

        return result;

    }

    /**
     * Gets the observed nucleotide for the given allele.
     * @param allele The allele whose nucleotide you are looking for.
     * @return The nucleotide, a character from the set {'A','C','G','T'}.
     */
    abstract public char getObservationForAllele(String allele);

    public abstract Collection<String> getAlleles();
}
