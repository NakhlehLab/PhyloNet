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

package edu.rice.cs.bioinfo.programs.phmm.src.optimize;

/**
 * Directly handle model manipulation in here.
 * Lightweight. Fiddles with SwitchingFrequencyRatio as a simple positive rate (member of a ratio).
 */

import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm.SwitchingFrequencyRatioTerm;
import edu.rice.cs.bioinfo.programs.phmm.src.runHmm.runHmm;

public class SwitchingFrequencyRatioTermParameter extends Parameter {
    // allow this to get basically to zero
    // divide the max by 10^10
    //public static final double DEFAULT_MINIMUM_RATE = MultivariateOptimizer.ABSOLUTE_ACCURACY;
    // careful - max depends on SwitchingFrequencyRatioTerm.getNumAlternatives()
    // eh - just go an order of magnitude below max
    //public static final double DEFAULT_INITIAL_RATE = MultivariateOptimizer.DEFAULT_INITIAL_RATE;
    // Like in original algorithm,
    // switching frequency approaching half or more can cause a degenerate corner case.
    // Lots of switching artifacts.
    //
    // Cap this so that maximum switching frequency is 0.25 .
    //
    // For two parental tree case, term is maximum  1/3 .
    //
    // Without any complicated structure to the ratio, the maximum is 1/(3(k-1)) where k is 
    // the number of alternative "states" (number of parental trees, or gene genealogies, or whatever).
    //
    // Conservatively set any non-self-transition probability to be at most 1/(3(k-1)),
    // ensuring that the self-transition probability is at least 3/4.
    //public static final double DEFAULT_MAXIMUM_RATE = MultivariateOptimizer.DEFAULT_MAXIMUM_RATE;

    // corresponds to maximum total non-self-transition switching frequency of 1/3
    //
    // actual constraint is applied in runHmm.calculateSwitchingFrequencies(...)
    //public static final double DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_TERMS_TOTAL_WEIGHT = 1.0 / 2.0;

    // rest calculated off of maximum
    //public static final double DEFAULT_MINIMUM_RATIO_TERM = 1e-6;
    public static final double DEFAULT_INITIAL_RATIO_TERM = 1.0;
    // any more than this will certainly cause normalization+rescaling in runHmm.calculateSwitchingFrequencies(...)
    //public static final double DEFAULT_MAXIMUM_RATIO_TERM = DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_TERMS_TOTAL_WEIGHT;
    //public static final double DEFAULT_MAXIMUM_RATIO_TERM = MultivariateOptimizer.DEFAULT_MAXIMUM_RATE;

    protected runHmm runHmmObject; // for pushing updated transition probabilities
    protected SwitchingFrequencyRatioTerm sfrt;

    protected double minimumRatioTerm;
    protected double maximumRatioTerm;

    /**
     * User must explicitly request whether or not an update is requested 
     * during construction.
     */
    public SwitchingFrequencyRatioTermParameter (String inName, 
						 double inValue, 
						 runHmm inRunHmmObject,
						 SwitchingFrequencyRatioTerm inSfrt,
						 boolean checkValueMinimumConstraintFlag,
						 boolean checkValueMaximumConstraintFlag,
						 boolean updateModelStateFlag,
						 double inMinimumRatioTerm,
						 double inMaximumRatioTerm
						 ) {
	// order forced by Java language constraints
	// need to delay the setValue() check until after runHmmObject reference ready
	super(inName, inValue, checkValueMinimumConstraintFlag, checkValueMinimumConstraintFlag, false);

	this.runHmmObject = inRunHmmObject;
	this.sfrt = inSfrt;

	if (updateModelStateFlag) {
	    updateModelState();
	}

        setMinimumValue(inMinimumRatioTerm);
	setMaximumValue(inMaximumRatioTerm);
    }

    // really, can go to zero
    public double getMinimumValue () {
	return (minimumRatioTerm);
    }

    public double getDefaultInitialValue () {
	return (DEFAULT_INITIAL_RATIO_TERM);
    }

    public double getMaximumValue () {
	return (maximumRatioTerm);
    }

    public void setMinimumValue (double inMinimumRatioTerm) {
	minimumRatioTerm = inMinimumRatioTerm;
    }

    public void setMaximumValue (double inMaximumRatioTerm) {
	maximumRatioTerm = inMaximumRatioTerm;
    }

    public void updateModelState () {
	// invalidate caches
	// caches re-populated in updateTransitionProbabilities call below
	//
	// push to SwitchingFrequencyRatioTerm.setValue()
	sfrt.setValue(getValue());

	// propagate hybridization/etc. frequency parameters on to
	// HMM transition probability matrix
	//
	// Throw in caching into runHmm?
	// If parental branch lengths change, no need to 
	// re-calculate switching frequencies from switching-frequency-ratio-terms.
	runHmmObject.updateTransitionProbabilities();
    }

}
