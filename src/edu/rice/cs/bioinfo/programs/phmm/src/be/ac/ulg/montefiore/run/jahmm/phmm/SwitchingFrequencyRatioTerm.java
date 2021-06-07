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
 * Simple container class for ratio term (positive real number)
 * corresponding to a switching frequency parameter used to calculate
 * hidden state transition probabilities.
 *
 * No storage of which trees (parental tree or gene genealogy)
 * this corresponds to. Maintained externally using a BidirectionalMultimap.
 *
 * See writeup for details.
 */

package edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm;

import edu.rice.cs.bioinfo.programs.phmm.src.optimize.CalculationCache;

public class SwitchingFrequencyRatioTerm {
    protected String name;
    protected double value;
    // perform cache invalidation in here
    protected CalculationCache calculationCache;
    // set to true to invalidate cacheParentalTreeSwitchingFrequencyMap
    // otherwise invalidates cacheGeneGenealogySwitchingFrequencyMap
    //protected boolean invalidateParentalTreeSwitchingFrequencyMapFlag;
    // bleh, also need to store the maximum number of alternatives
    // (parental trees, gene genealogies, or whatever)
    // to help cap the switching frequencies
    protected int numAlternatives;
    
    // boolean inInvalidateParentalTreeSwitchingFrequencyMapFlag, 
    
    public SwitchingFrequencyRatioTerm (String inName, double inValue, CalculationCache inCalculationCache, int inNumAlternatives) {
	this.calculationCache = inCalculationCache;
	//this.invalidateParentalTreeSwitchingFrequencyMapFlag = inInvalidateParentalTreeSwitchingFrequencyMapFlag;
	this.numAlternatives = inNumAlternatives;

	setName(inName);
	setValue(inValue);
    }
    
    public String getName () {
	return (name);
    }
    
    public double getValue () {
	return (value);
    }

    public void setName (String inName) {
	name = inName;
    }

    public void setValue (double inValue) {
	if (inValue > 0.0) {
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
	    throw (new RuntimeException("ERROR: SwitchingFrequencyRatioTerm.setValue(double) called with non-positive number."));
	}
    }

    public int getNumAlternatives () {
	return (numAlternatives);
    }

    /**
     * For map support.
     */
    public boolean equals (Object obj) {
	if (!(obj instanceof SwitchingFrequencyRatioTerm)) {
	    return (false);
	}

	SwitchingFrequencyRatioTerm so = (SwitchingFrequencyRatioTerm) obj;
	return (this.getName().equals(so.getName()));
    }

    public int hashCode() {
	return (this.getName().hashCode());
    }

    public String toString () {
	return (getName() + ": " + getValue());
    }
}
