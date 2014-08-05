package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations;

import java.util.*;

public class ManyNucleotideObservation extends NucleotideObservation
{


    Map<String,String> dna;
    int index;

    int[] nucleotides;

    public ManyNucleotideObservation(Map<String,String> dna, int index)
    {
        this.dna = dna;
        this.index = index;


        nucleotides = new int[dna.entrySet().size()];
        int i = 0;
        for (String value : dna.values())
        {
            nucleotides[i++] = value.charAt(index);
        }
    }

    @Override
    public char getObservationForAllele(String allele) {
        return dna.get(allele).charAt(index);
    }

    @Override
    public Collection<String> getAlleles() {
        return dna.keySet();
    }


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ManyNucleotideObservation that = (ManyNucleotideObservation) o;

        return Arrays.equals(nucleotides,that.nucleotides);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(nucleotides);
    }

    @Override
    public String toString() {
        String dnaString = "{";
        for (String key : dna.keySet())
        {
            dnaString += key + "=" + dna.get(key).charAt(index) + ",";
        }

        return "ManyNucleotideObservation{" +
                "dna=" + dnaString +
                '}';
    }

    public static void main(String[] args)
    {
        Map<String,String> dnaAllel = new HashMap<String,String>();
        dnaAllel.put("A","ABCD");
        dnaAllel.put("B","AABB");
        dnaAllel.put("C","ABAB");

        List<NucleotideObservation> obs = NucleotideObservation.getObservations(dnaAllel);

        for (NucleotideObservation ob : obs)
        {
            System.out.println(ob + " " + Arrays.toString(((ManyNucleotideObservation) ob).nucleotides) + " " + ob.hashCode());
        }
    }
}
