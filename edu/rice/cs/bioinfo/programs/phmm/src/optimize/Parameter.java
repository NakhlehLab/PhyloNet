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
 * Simple state container for length parameter.
 * Reduce indirection between IDs and objects a bit.
 *
 * Don't clone.
 */

// make a constrained subclass that 
// takes in reference to a set of length parameters
// and a ParameterConstraintSet

package edu.rice.cs.bioinfo.programs.phmm.src.optimize;

import java.util.Stack;

public abstract class Parameter extends ParameterConstraintSet {
    protected Stack<Double> cacheValues;

    public Parameter (String inName) {
	super(inName);
	cacheValues = new Stack<Double>();
    }
    
    // public Parameter (String inName, double inValue) {
    // 	super(inName, inValue);
    // }

    // force user to think carefully about this
    // alternative constructor with relaxed constraints
    public Parameter (String inName, 
		      double inValue,
		      boolean checkValueMinimumConstraintFlag,
		      boolean checkValueMaximumConstraintFlag,
		      boolean updateModelStateFlag) {
	super(inName);
	cacheValues = new Stack<Double>();
	setValue(inValue, checkValueMinimumConstraintFlag, checkValueMaximumConstraintFlag, updateModelStateFlag);
    }

    public void setName (String inName) {
	super.setName(inName);
    }

    /**
     * WARNING - will push updated state to model!
     */
    public void setValue (double inValue) {
	setValue(inValue, true, true, true);
    }

    /**
     * WARNING - don't set updateModelStateFlag to false unless you know what you're doing.
     *
     * Parameter object's value is now coupled to 
     * HMM state only by updateModelState(), as controlled
     * by updateModelStateFlag.
     */
    public void setValue (double inValue, 
			  boolean checkValueMinimumConstraintFlag,
			  boolean checkValueMaximumConstraintFlag,
			  boolean updateModelStateFlag) {
	// paranoid
	if (checkValueMinimumConstraintFlag && (inValue < getMinimumValue())) {
	    throw (new RuntimeException("ERROR: Parameter.setValue(...) called with value smaller than getMinimumValue() amount. " + getName() + " " + inValue + " " + getMinimumValue()));
	}

	// kliu - this can cause an issue for FrequencyParameter
	// constructor
	if (checkValueMaximumConstraintFlag && (inValue > getMaximumValue())) {
	    throw (new RuntimeException("ERROR: Parameter.setValue(...) called with value greater than getMaximumValue() amount. " + getName() + " " + inValue + " " + getMaximumValue()));
	}

	super.setValue(inValue);

	if (updateModelStateFlag) {
	    updateModelState();
	}
    }

    /**
     * Also provide cache and restore feature.
     */
    public void pushCacheValue () {
	cacheValues.push(new Double(getValue()));
    }
    
    /**
     * WARNING - don't call prior to calling cacheValue().
     * 
     * Add more flags for additional control. See setValue(...).
     * Set setValueFlag to false to only pop, no setValue(...) call.
     */
    public void popCacheValue (boolean checkValueMinimumConstraintFlag,
			       boolean checkValueMaximumConstraintFlag,
			       boolean updateModelStateFlag,
			       boolean setValueFlag) {
	if (setValueFlag) {
	    // This will call updateModelState() appropriately.
	    setValue(cacheValues.pop().doubleValue(),
		     checkValueMinimumConstraintFlag,
		     checkValueMaximumConstraintFlag,
		     updateModelStateFlag);
	}
	else {
	    cacheValues.pop();
	}
    }

    /**
     * WARNING - don't call prior to calling cacheValue().
     */
    // public void popCacheValue () {
    // 	popCacheValue(true);
    // }

    // kliu - also need min/max
    public abstract double getMinimumValue ();
    public abstract double getDefaultInitialValue ();
    public abstract double getMaximumValue ();

    // Migrate model updates to Parameter class and subclasses.
    public abstract void updateModelState ();
    
    // derp
    public boolean equals (Object obj) {
	if (!(obj instanceof Parameter)) {
	    return (false);
	}
	Parameter lp = (Parameter) obj;
	return (this.getName().equals(lp.getName()));
    }
}
