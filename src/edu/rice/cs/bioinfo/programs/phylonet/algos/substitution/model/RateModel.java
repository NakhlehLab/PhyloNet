package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model;

import jeigen.DenseMatrix;

/**
 * An interface for nucleotide rate models, with certain rates of transitions from one base to another.
 */
public interface RateModel
{
    /**
     * Returns the equilibrium column vector for this rate model.
     * The following relationship should hold: rateMatrix * equalibriumVector = 0.
     * @return The equilibrium positions as a column vector.
     */
    DenseMatrix getEquilibriumVector();

    /**
     * Returns the rate matrix for this model.
     * The sum of each column should be 0.
     * @return The rate matrix.
     */
    DenseMatrix getRateMatrix();
}
