package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.substitution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.substitution.EigenDecomposition;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.substitution.Frequencies;

/**
 * Created by wendingqiao on 5/3/16.
 * Specifies substitution model from which a transition probability matrix for a given distance can be obtained.
 */
public interface SubstitutionModel {

    /**
     * Gets the transition probability matrix where distance = (startTime - endTime) * rate.
     */
    void getTransitionProbabilities(double startTime, double endTime, double rate, double[] matrix);

    /**
     * Gets matrix Q.
     */
    double[][] getRateMatrix();

    double[] getFrequencies();

    int getStateCount();


    /**
     * Gets the Eigen decomposition of matrix Q.
     */
    EigenDecomposition getEigenDecomposition();

    /**
     * returns whether substitution model can return complex diagonalizations.
     */
    boolean canReturnComplexDiagonalization();

    /**
     * basic implementation of substitution model
     */
    abstract class Base extends StateNode implements SubstitutionModel {

        protected Frequencies _frequencies;
        protected int _nrOfStates;

        public Base(Frequencies freq) {
            _frequencies = freq;
        }

        @Override
        public double[] getFrequencies() {
            return _frequencies.getFreqs();
        }

        @Override
        public int getStateCount() {
            return _nrOfStates;
        }

        @Override
        public double[][] getRateMatrix() {
            return null;
        }

        @Override
        public boolean canReturnComplexDiagonalization() {
            return false;
        }

    }
}
