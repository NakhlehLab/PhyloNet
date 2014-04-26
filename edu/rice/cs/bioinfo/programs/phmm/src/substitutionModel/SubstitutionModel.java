/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
