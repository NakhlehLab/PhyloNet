package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model;

import be.ac.ulg.montefiore.run.jahmm.Observation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * An observation class which represents seeing one nucleotide per allele.
 */
public class JahmmNucleotideObservation extends Observation
{

    NucleotideObservation self;
    public JahmmNucleotideObservation(NucleotideObservation self)
    {
        this.self = self;
    }

    public static List<JahmmNucleotideObservation> wrapObservations(List<NucleotideObservation> nucObs)
    {
        List<JahmmNucleotideObservation> result = new ArrayList<JahmmNucleotideObservation>();
        for (NucleotideObservation obs: nucObs)
            result.add(new JahmmNucleotideObservation(obs));
        return result;
    }

    public NucleotideObservation getObservation()
    {
        return self;
    }

    @Override
    public String toString(NumberFormat numberFormat)
    {
        return toString();
    }

    @Override
    public String toString()
    {
        return self.toString();
    }
}
