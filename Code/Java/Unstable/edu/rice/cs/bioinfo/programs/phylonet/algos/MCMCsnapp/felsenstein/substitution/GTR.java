package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.substitution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by wendingqiao on 5/3/16.
 */
public class GTR extends GeneralSubstitutionModel {

    // rates should be set up beforehand
    // default: all 1
    public GTR(Frequencies freq, double[] rates) {
        super(freq, new double[] {
                rates[0], // A->C
                rates[1], // A->G
                rates[2], // A->T
                rates[0], // C->A
                rates[3], // C->G
                rates[4], // C->T
                rates[1], // G->A
                rates[3], // G->C
                rates[5], // G->T
                rates[2], // T->A
                rates[4], // T->C
                rates[5], // T->G
        });
    }

    public GTR(Frequencies freq) {
        this(freq, new double[]{1,1,1,1,1,1,1,1,1,1,1,1});
    }

    @Override
    public double propose() {
        return 0;
    }

    @Override
    public void undo() {

    }

    @Override
    public void accept() {

    }

    @Override
    public void reject() {

    }

    @Override
    public double logDensity() {
        return 0;
    }

    @Override
    public boolean mayViolate() {
        return false;
    }

    @Override
    public boolean isValid() {
        return true;
    }

    // test
    public static void main(String[] args) {
        Map<String, String> map = new HashMap<>();
        map.put("A", "ATCG");
        map.put("B", "ATTG");
        map.put("C", "AGAG");
        MarkerSeq aln = new MarkerSeq(map);
        Frequencies freq = new Frequencies(aln, false);
        double[] matrix = new double[16];
        GTR gtr = new GTR(freq);
        gtr.getTransitionProbabilities(1.0, 0.0, 0.1, matrix);
        System.out.println(Arrays.toString(gtr.getFrequencies()));
        System.out.println(Arrays.toString(matrix));
    }
}
