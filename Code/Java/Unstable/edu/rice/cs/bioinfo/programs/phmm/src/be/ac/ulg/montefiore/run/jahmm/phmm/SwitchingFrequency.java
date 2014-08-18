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
 * Simple container class for switching frequency.
 *
 * See writeup for details.
 */

package edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm;

import edu.rice.cs.bioinfo.programs.phmm.src.optimize.CalculationCache;

public class SwitchingFrequency extends SwitchingFrequencyRatioTerm {
    public static final String GAMMA = "gamma";

    public SwitchingFrequency (String inName, double inValue, CalculationCache inCalculationCache, int inNumAlternatives) {
	super(inName, inValue, inCalculationCache, inNumAlternatives);
    }

    public void setValue (double inValue) {
	if ((inValue >= 0.0) & (inValue <= 1.0)) {
	    value = inValue;
	    // wipe cache
	    //calculationCache.cacheSwitchingFrequencyMap = null;

	    // if (invalidateParentalTreeSwitchingFrequencyMapFlag) {
	    // 	calculationCache.cacheParentalTreeSwitchingFrequencyMap = null;
	    // }
	    // else {
	    // 	calculationCache.cacheGeneGenealogySwitchingFrequencyMap = null;
	    // }
	}
	else {
	    throw (new RuntimeException("ERROR: SwitchingFrequency.setValue(double) called with argument outside of [0,1] range."));
	}
    }

}

