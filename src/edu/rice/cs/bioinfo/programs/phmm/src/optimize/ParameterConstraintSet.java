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
 * Just need some storage for constraint.
 * Immutable after creation.
 *
 * Use a BidirectionalMultimap to assign Parameter objects to a ParameterConstraintSet
 */

package edu.rice.cs.bioinfo.programs.phmm.src.optimize;

public class ParameterConstraintSet {
    protected String name;
    // kliu - don't get too fancy with generic type for value
    // Apache's BrentOptimizer implementation only permits a double for parameter values.
    protected double value;

    public ParameterConstraintSet (String inName) {
    	setName(inName);
    }

    public ParameterConstraintSet (String inName, double inValue) {
	setName(inName);
	setValue(inValue);
    }

    public String getName () {
	return (name);
    }

    public double getValue () {
	return (value);
    }

    // immutable once created
    protected void setName (String inName) {
	name = inName;
    }

    protected void setValue (double inValue) {
	value = inValue;
    }

    // ParameterConstraintSet objects are hashCode() and equals() equivalent to
    // their name Strings

    public String toString () {
	return (name.toString() + ": " + Double.toString(value));
    }

    public int hashCode () {
	return (name.hashCode());
    }

    // derp
    public boolean equals (Object obj) {
	if (!(obj instanceof ParameterConstraintSet)) {
	    return (false);
	}
	ParameterConstraintSet lpcs = (ParameterConstraintSet) obj;
	return (this.getName().equals(lpcs.getName()));
    }
}

