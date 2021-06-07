package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.StateNode;

/**
 * Created by wendingqiao on 5/3/16.
 * Specifies transition probability matrix with no restrictions on the rates other
 * than that one of the is equal to one and the others are specified relative to
 * this unit rate. Works for any number of states
 */
public class GeneralSubstitutionModel extends SubstitutionModel.Base {

    protected double[] _transitionRates; // AC, AG, AT, CG, CT, GT

    protected double[][] _rateMatrix;

    protected double[] _relativeRates;
    protected double[] _storedRelativeRates;

    protected EigenSystem _eigenSystem;

    protected EigenDecomposition _eigenDecomposition;
    private EigenDecomposition _storedEigenDecomposition;

    protected boolean _updateMatrix = true;
    private boolean _storedUpdateMatrix = true;

    public GeneralSubstitutionModel(Frequencies freq, double[] transitionRates) {
        super (freq);
        _transitionRates = transitionRates;
        _updateMatrix = true;
        _nrOfStates = _frequencies.getFreqs().length;
        if (_transitionRates.length != _nrOfStates * (_nrOfStates - 1)) {
            throw new IllegalArgumentException("Dimension of transition rates is " + _transitionRates.length);
        }
        _eigenSystem = new EigenSystem(_nrOfStates);
        _rateMatrix = new double[_nrOfStates][_nrOfStates];
        _relativeRates = new double[_transitionRates.length];
        _storedRelativeRates = new double[_transitionRates.length];
    }

    @Override
    public void getTransitionProbabilities(double startTime, double endTime, double rate, double[] matrix) {
        double distance = (startTime - endTime) * rate;

        int i, j, k;
        double temp;

        synchronized (this) {
            if (_updateMatrix) {
                setupRelativeRates();
                setupRateMatrix();
                _eigenDecomposition = _eigenSystem.decomposeMatrix(_rateMatrix);
                _updateMatrix = false;
            }
        }

        double[] iexp = new double[_nrOfStates * _nrOfStates];
        double[] Evec = _eigenDecomposition.getEigenVectors();
        double[] Ievc = _eigenDecomposition.getInverseEigenVectors();
        double[] Eval = _eigenDecomposition.getEigenValues();
        for (i = 0; i < _nrOfStates; i++) {
            temp = Math.exp(distance * Eval[i]);
            for (j = 0; j < _nrOfStates; j++) {
                iexp[i * _nrOfStates + j] = Ievc[i * _nrOfStates + j] * temp;
            }
        }

        int u = 0;
        for (i = 0; i < _nrOfStates; i++) {
            for (j = 0; j < _nrOfStates; j++) {
                temp = 0.0;
                for (k = 0; k < _nrOfStates; k++) {
                    temp += Evec[i * _nrOfStates + k] * iexp[k * _nrOfStates + j];
                }

                matrix[u] = Math.abs(temp);
                u++;
            }
        }
    }

    public double[][] getRateMatrix() {
        return _rateMatrix.clone();
    }

    protected void setupRelativeRates() {
        _relativeRates = _transitionRates.clone();
    }

    protected void setupRateMatrix() {
        double[] freqs = _frequencies.getFreqs();
        for (int i = 0; i < _nrOfStates; i++) {
            _rateMatrix[i][i] = 0;
            for (int j = 0; j < i; j++) {
                _rateMatrix[i][j] = _relativeRates[i * (_nrOfStates - 1) + j];
            }
            for (int j = i + 1; j < _nrOfStates; j++) {
                _rateMatrix[i][j] = _relativeRates[i * (_nrOfStates - 1) + j - 1];
            }
        }
        // bring in frequencies
        for (int i = 0; i < _nrOfStates; i++) {
            for (int j = i + 1; j < _nrOfStates; j++) {
                _rateMatrix[i][j] *= freqs[j];
                _rateMatrix[j][i] *= freqs[i];
            }
        }
        // set up diagonal
        for (int i = 0; i < _nrOfStates; i++) {
            double sum = 0.0;
            for (int j = 0; j < _nrOfStates; j++) {
                if (i != j)
                    sum += _rateMatrix[i][j];
            }
            _rateMatrix[i][i] = -sum;
        }
        // normalize rate matrix to one expected substitution per unit time
        double subst = 0.0;
        for (int i = 0; i < _nrOfStates; i++)
            subst += -_rateMatrix[i][i] * freqs[i];

        for (int i = 0; i < _nrOfStates; i++) {
            for (int j = 0; j < _nrOfStates; j++) {
                _rateMatrix[i][j] = _rateMatrix[i][j] / subst;
            }
        }
    }

    @Override
    public EigenDecomposition getEigenDecomposition() {
        synchronized (this) {
            if (_updateMatrix) {
                setupRelativeRates();
                setupRateMatrix();
                _eigenDecomposition = _eigenSystem.decomposeMatrix(_rateMatrix);
                _updateMatrix = false;
            }
        }
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
}
