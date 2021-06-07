package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model;

import be.ac.ulg.montefiore.run.jahmm.Opdf;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.NucleotideProbabilityAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.text.NumberFormat;
import java.util.Collection;

public class HmmState implements Opdf<JahmmNucleotideObservation>{

    Tree speciesTree;
    Tree geneTree;
    SubstitutionModel model;
    NucleotideProbabilityAlgorithm nucleotideProbabilityAlgorithm;

    public HmmState(Tree speciesTree, Tree geneTree,SubstitutionModel model, NucleotideProbabilityAlgorithm nucleotideProbabilityAlgorithm){
        this.speciesTree = speciesTree;
        this.geneTree = geneTree;
        this.model = model;
        this.nucleotideProbabilityAlgorithm = nucleotideProbabilityAlgorithm;
    }


    @Override
    public double probability(JahmmNucleotideObservation nucleotideObservation) {
        return nucleotideProbabilityAlgorithm.getProbabilityCached(nucleotideObservation.getObservation());
    }

    @Override
    public JahmmNucleotideObservation generate() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void fit(JahmmNucleotideObservation... oa) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void fit(JahmmNucleotideObservation[] o, double[] weights) {
        throw new UnsupportedOperationException();
    }


    @Override
    public void fit(Collection co) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void fit(Collection co, double[] weights) {
        throw new UnsupportedOperationException();
    }

    @Override
    public String toString(NumberFormat numberFormat) {
        return toString();
    }



    @Override
    public Opdf<JahmmNucleotideObservation> clone() {
        throw new UnsupportedOperationException();

    }

    @Override
    public String toString() {
        return "HmmState{" +
                "speciesTree=" + speciesTree +
                ", geneTree=" + geneTree +
                '}';
    }
}
