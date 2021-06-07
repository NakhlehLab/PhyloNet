package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations;

import java.text.NumberFormat;
import java.util.Collection;
import java.util.Map;

public class OneNucleotideObservation extends NucleotideObservation
{

    Map<String,Character> nucleotides;

    public OneNucleotideObservation(Map<String,Character> nucleotides)
    {
        this.nucleotides = nucleotides;
    }

    public char getObservationForAllele(String allele)
    {
        return nucleotides.get(allele);
    }

    @Override
    public Collection<String> getAlleles() {
        return nucleotides.keySet();
    }

    @Override
    public String toString() {
        return "OneNucleotideObservation{" +
                "nucleotides=" + nucleotides +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        OneNucleotideObservation that = (OneNucleotideObservation) o;

        if (!nucleotides.equals(that.nucleotides)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return nucleotides.hashCode();
    }



}
