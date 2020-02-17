package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.substitution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.datatype.DataType;

import java.util.Arrays;

/**
 * Created by wendingqiao on 5/3/16.
 */
public class Frequencies extends StateNode {

    private static final double MIN_VAL = 1.0E-10;

    protected Alignment _sequences; // Sequence data for which frequencies are calculated
    protected boolean _estimate = false; // Whether to estimate the frequencies from data or not
    protected double[] _frequencies; // A set of frequencies specified as space separated values summing to 1
    boolean _needsUpdate;

    public Frequencies(Alignment aln, boolean estimate) {
        _sequences = aln;
        _estimate = estimate;
        // initialization - equally distributed
        int states = _sequences.getMaxStateCount();
        _frequencies = new double[states];
        for (int i = 0; i < states; i++) {
            _frequencies[i] = 1.0 / states;
        }
        if (_estimate) {
            estimateFrequencies();
            checkFrequencies();
        }
    }

    public Frequencies(Alignment aln, double[] frequencies) {
        double sum = getSumOfFrequencies(frequencies);
        if (Math.abs(sum - 1.0) > 1e-6) {
            throw new IllegalArgumentException("Frequencies do not add up to 1");
        }
        _sequences = aln;
        _frequencies = frequencies;
    }

    /**
     * return up to date frequencies
     */
    public double[] getFreqs() {
        return _frequencies.clone();
    }

    /**
     * Estimate from sequence alignment.
     * This version matches the implementation in Beast 1 & PAUP  *
     */
    void estimateFrequencies() {
        DataType dataType = _sequences.getDataType();
        int stateCount = _sequences.getMaxStateCount();

        int attempts = 0;
        double difference;
        do {
            double[] tmpFreq = new double[stateCount];

            double total = 0.0;
            for (int i = 0; i < _sequences.getPatternCount(); i++) {
                int[] pattern = _sequences.getPattern(i);
                double weight = _sequences.getPatternWeight(i);

                for (int value : pattern) {
                    int[] codes = dataType.getStatesForCode(value);
                    double sum = 0.0;
                    for (int codeIndex : codes) {
                        sum += _frequencies[codeIndex];
                    }
                    for (int codeIndex : codes) {
                        double tmp = (_frequencies[codeIndex] * weight) / sum;
                        tmpFreq[codeIndex] += tmp;
                        total += tmp;
                    }
                }
            }
            difference = 0.0;
            for (int i = 0; i < stateCount; i++) {
                difference += Math.abs((tmpFreq[i] / total) - _frequencies[i]);
                _frequencies[i] = tmpFreq[i] / total;
            }
            attempts++;
        } while (difference > 1E-8 && attempts < 1000);
        System.out.println("Starting frequencies: " + Arrays.toString(_frequencies));
    }

    /**
     * Ensures that frequencies are not smaller than MINFREQ and
     * that two frequencies differ by at least 2*MINFDIFF.
     * This avoids potential problems later when eigenvalues
     * are computed.
     */
    private void checkFrequencies() {
        int maxi = 0;
        double sum = 0.0;
        double maxfreq = 0.0;
        for (int i = 0; i < _frequencies.length; i++) {
            double freq = _frequencies[i];
            if (freq < MIN_VAL) _frequencies[i] = MIN_VAL;
            if (freq > maxfreq) {
                maxfreq = freq;
                maxi = i;
            }
            sum += _frequencies[i];
        }
        double diff = 1.0 - sum;
        _frequencies[maxi] += diff;

        for (int i = 0; i < _frequencies.length - 1; i++) {
            for (int j = i + 1; j < _frequencies.length; j++) {
                if (_frequencies[i] == _frequencies[j]) {
                    _frequencies[i] += MIN_VAL;
                    _frequencies[j] -= MIN_VAL;
                }
            }
        }
    }

    /**
     * @param frequencies the frequencies
     * @return return the sum of frequencies
     */
    private double getSumOfFrequencies(double[] frequencies) {
        double total = 0.0;
        for (double frequency : frequencies) {
            total += frequency;
        }
        return total;
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
}