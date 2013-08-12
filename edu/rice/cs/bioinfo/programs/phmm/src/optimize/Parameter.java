/**
 * Simple state container for length parameter.
 * Reduce indirection between IDs and objects a bit.
 *
 * Don't clone.
 */

// make a constrained subclass that 
// takes in reference to a set of length parameters
// and a ParameterConstraintSet

package optimize;

public abstract class Parameter extends ParameterConstraintSet {
    protected double cachedValue;

    public Parameter (String inName, double inValue) {
	super(inName, inValue);
    }

    // alternative constructor with relaxed constraints
    public Parameter (String inName, 
		      double inValue,
		      boolean checkValueMinimumConstraintFlag,
		      boolean checkValueMaximumConstraintFlag) {
	super(inName);
	setValue(inValue, checkValueMinimumConstraintFlag, checkValueMaximumConstraintFlag);
    }

    public void setName (String inName) {
	super.setName(inName);
    }

    public void setValue (double inValue) {
	setValue(inValue, true, true);
    }

    public void setValue (double inValue, 
			  boolean checkValueMinimumConstraintFlag,
			  boolean checkValueMaximumConstraintFlag) {
	// paranoid
	if (checkValueMinimumConstraintFlag && (inValue < getMinimumValue())) {
	    throw (new RuntimeException("ERROR: Parameter.setValue(...) called with value smaller than getMinimumValue() amount."));
	}

	// kliu - this can cause an issue for FrequencyParameter
	// constructor
	if (checkValueMaximumConstraintFlag && (inValue > getMaximumValue())) {
	    throw (new RuntimeException("ERROR: Parameter.setValue(...) called with value greater than getMaximumValue() amount."));
	}

	super.setValue(inValue);
    }

    /**
     * Also provide cache and restore feature.
     */
    public void cacheValue () {
	cachedValue = getValue();
    }
    
    /**
     * WARNING - don't call prior to calling cacheValue().
     */
    public void restoreCachedValue () {
	setValue(cachedValue);
    }
    
    // kliu - also need min/max
    public abstract double getMinimumValue ();
    public abstract double getDefaultInitialValue ();
    public abstract double getMaximumValue ();
    
    // derp
    public boolean equals (Object obj) {
	if (!(obj instanceof Parameter)) {
	    return (false);
	}
	Parameter lp = (Parameter) obj;
	return (this.getName().equals(lp.getName()));
    }
}
