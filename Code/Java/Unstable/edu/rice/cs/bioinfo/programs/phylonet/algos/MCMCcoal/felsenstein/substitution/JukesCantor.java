package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.substitution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.substitution.EigenDecomposition;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.substitution.Frequencies;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.substitution.SubstitutionModel;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by wendingqiao on 5/3/16.
 * Jukes Cantor substitution model: all rates equal and uniformly distributed frequencies
 */
public class JukesCantor extends SubstitutionModel.Base {

    EigenDecomposition _eigenDecomposition;

    public JukesCantor(Frequencies freq) {
        super(freq);
        double[] eval = new double[]{0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333};
        double[] evec = new double[]{1.0, 2.0, 0.0, 0.5, 1.0, -2.0, 0.5, 0.0, 1.0, 2.0, 0.0, -0.5, 1.0, -2.0, -0.5, 0.0};
        double[] ivec = new double[]{0.25, 0.25, 0.25, 0.25, 0.125, -0.125, 0.125, -0.125, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0, 0.0};
        _eigenDecomposition = new EigenDecomposition(evec, ivec, eval);

        if (freq._estimate) {
            throw new RuntimeException("Frequencies must not be specified in Jukes-Cantor model. They are assumed equal.");
        }
    }

    @Override
    public double[] getFrequencies() {
        return super.getFrequencies();
    }

    @Override
    public void getTransitionProbabilities(double startTime, double endTime, double rate, double[] matrix) {
        double delta = 4.0 / 3.0 * (startTime - endTime);
        double pStay = (1.0 + 3.0 * Math.exp(-delta * rate)) / 4.0;
        double pMove = (1.0 - Math.exp(-delta * rate)) / 4.0;
        // fill the matrix with move probabilities
        Arrays.fill(matrix, pMove);
        // fill the diagonal
        for (int i = 0; i < 4; i++) {
            matrix[i * 5] = pStay;
        }
    }

    @Override
    public EigenDecomposition getEigenDecomposition() {
        return _eigenDecomposition;
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
        Alignment aln = new Alignment(map);
        Frequencies freq = new Frequencies(aln, false);
        JukesCantor jc = new JukesCantor(freq);
        double[] matrix = new double[16];
        jc.getTransitionProbabilities(1.0, 0.0, 0.1, matrix);
        System.out.println(Arrays.toString(jc.getFrequencies()));
        System.out.println(Arrays.toString(matrix));
    }
}