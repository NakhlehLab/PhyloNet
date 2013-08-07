/**
 * Calculate a substitution probability matrix
 * from a sequence evolution model, e.g. GTR.
 */

package substitutionModel;

import util.*;

public interface SubstitutionModel {
    /**
     * To support write once on construction, and
     * then read forever after. 
     * Protect internal state from external mutation.
     * Since annoying to propagate mutated state to external dependent
     * calculated quantities.
     */
    public SubstitutionModel deepCopyClone ();

    public Alphabet getAlphabet ();

    /**
     * Check for reversibility.
     * Some models are that way by construction,
     * others will have to check their parameters to see
     * if reversibility applies for that set of parameter choices.
     */
    public boolean checkReversible ();

    /**
     * From Felsenstein 2001 Inferring Phylogenies
     * R[t] = e^{Qt}
     * R the substitution probability matrix
     * Q the substitution rate matrix
     * t the branch length
     */
    public double[][] calculateProbabilitiesFromRates (double time);

    /**
     * Also need stationary probabilities.
     */
    public double[] getStationaryProbabilities ();
}
