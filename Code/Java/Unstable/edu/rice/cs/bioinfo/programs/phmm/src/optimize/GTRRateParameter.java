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

import substitutionModel.GTRSubstitutionModel;

public class GTRRateParameter extends Parameter {
    public static final double DEFAULT_MINIMUM_RATE = MultivariateOptimizer.DEFAULT_MINIMUM_RATE;
    public static final double DEFAULT_INITIAL_RATE = MultivariateOptimizer.DEFAULT_INITIAL_RATE;
    public static final double DEFAULT_MAXIMUM_RATE = MultivariateOptimizer.DEFAULT_MAXIMUM_RATE;

    protected GTRSubstitutionModel gtrSubstitutionModel;
    // choice of which parameter to optimize
    protected int index;
    protected CalculationCache calculationCache;

    public GTRRateParameter (String inName, 
			     double inValue,
			     GTRSubstitutionModel inGTRSubstitutionModel,
			     int inIndex,
			     CalculationCache inCalculationCache,
			     boolean checkValueMinimumConstraintFlag,
			     boolean checkValueMaximumConstraintFlag,
			     boolean updateModelStateFlag) {
	super(inName, inValue, checkValueMinimumConstraintFlag, checkValueMaximumConstraintFlag, false);

	this.gtrSubstitutionModel = inGTRSubstitutionModel;
	this.calculationCache = inCalculationCache;
	// paranoid
	if ((inIndex < 0) || (inIndex >= gtrSubstitutionModel.getRateParameterCount())) {
	    throw (new RuntimeException("ERROR: index in GTRRateParameter(...) is out of bounds. " + inIndex));
	}
	this.index = inIndex;
	
	if (updateModelStateFlag) {
	    updateModelState();
	}
    }

    public double getMinimumValue () {
	return (DEFAULT_MINIMUM_RATE);
    }

    public double getDefaultInitialValue () {
	return (DEFAULT_INITIAL_RATE);
    }

    public double getMaximumValue () {
	return (DEFAULT_MAXIMUM_RATE);
    }

    public void updateModelState () {
	double[] originalRateParameters = gtrSubstitutionModel.getOriginalRateParameters();
	// set by reference
	originalRateParameters[index] = this.getValue();
	gtrSubstitutionModel.updateRateMatrix();

	// totally clear out associated caches
	calculationCache.cacheSubstitutionProbabilityMatrix.clear();
	calculationCache.cacheSubstitutionProbability.clear();

    }
}
