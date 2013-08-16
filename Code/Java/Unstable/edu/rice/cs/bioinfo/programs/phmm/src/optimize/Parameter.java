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

import java.util.Stack;

public abstract class Parameter extends ParameterConstraintSet {
    protected Stack<Double> cacheValues;

    public Parameter (String inName, double inValue) {
    	super(inName, inValue);
    }

    // force user to think carefully about this
    // alternative constructor with relaxed constraints
    public Parameter (String inName, 
		      double inValue,
		      boolean checkValueMinimumConstraintFlag,
		      boolean checkValueMaximumConstraintFlag,
		      boolean updateModelStateFlag) {
	super(inName);
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
     */
    public void setValue (double inValue, 
			  boolean checkValueMinimumConstraintFlag,
			  boolean checkValueMaximumConstraintFlag,
			  boolean updateModelStateFlag) {
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

	if (updateModelStateFlag) {
	    updateModelState();
	}
    }

    /**
     * Also provide cache and restore feature.
     */
    public void pushCacheValue () {
	cacheValues.push(getValue());
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
	    setValue(cacheValues.pop(),
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
