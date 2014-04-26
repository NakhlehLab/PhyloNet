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

package optimize;

/**
 * Directly handle model manipulation in here.
 * Lightweight. Fiddles with SwitchingFrequency as a (capped) probability.
 */

import be.ac.ulg.montefiore.run.jahmm.phmm.SwitchingFrequency;
import runHmm.runHmm;

public class SwitchingFrequencyParameter extends Parameter {
    // Like in original algorithm,
    // switching frequency approaching half or more can cause a degenerate corner case.
    // Lots of switching artifacts.
    //
    // Cap this so that maximum switching frequency is 0.25 .
    //
    public static final double DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_FREQUENCY_TOTAL_WEIGHT = 0.25;

    // rest calculated off of maximum
    public static final double DEFAULT_MINIMUM_FREQUENCY = MultivariateOptimizer.DEFAULT_MINIMUM_PROBABILITY;
    public static final double DEFAULT_INITIAL_FREQUENCY = MultivariateOptimizer.DEFAULT_INITIAL_PROBABILITY;
    public static final double DEFAULT_MAXIMUM_FREQUENCY = DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_FREQUENCY_TOTAL_WEIGHT;

    protected runHmm runHmmObject; // for pushing updated transition probabilities
    protected SwitchingFrequency sf;

    /**
     * User must explicitly request whether or not an update is requested 
     * during construction.
     */
    public SwitchingFrequencyParameter (String inName, 
					double inValue, 
					runHmm inRunHmmObject,
					SwitchingFrequency inSf,
					boolean checkValueMinimumConstraintFlag,
					boolean checkValueMaximumConstraintFlag,
					boolean updateModelStateFlag) {
	// order forced by Java language constraints
	// need to delay the setValue() check until after runHmmObject reference ready
	super(inName, inValue, checkValueMinimumConstraintFlag, checkValueMinimumConstraintFlag, false);

	this.runHmmObject = inRunHmmObject;
	this.sf = inSf;

	if (updateModelStateFlag) {
	    updateModelState();
	}
    }

    // really, can go to zero
    public double getMinimumValue () {
	return (DEFAULT_MINIMUM_FREQUENCY);
    }

    public double getDefaultInitialValue () {
	return (DEFAULT_INITIAL_FREQUENCY);
    }

    public double getMaximumValue () {
	return (DEFAULT_MAXIMUM_FREQUENCY);
    }

    public void updateModelState () {
	// invalidate caches
	// caches re-populated in updateTransitionProbabilities call below
	//
	// push to SwitchingFrequency.setValue()
	sf.setValue(getValue());

	// propagate hybridization/etc. frequency parameters on to
	// HMM transition probability matrix
	runHmmObject.updateTransitionProbabilities();
    }

}
