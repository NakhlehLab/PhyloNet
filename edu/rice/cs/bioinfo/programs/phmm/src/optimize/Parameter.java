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
    public Parameter (String inName, double inValue) {
	super(inName, inValue);
    }

    public void setName (String inName) {
	super.setName(inName);
    }

    public void setValue (double inValue) {
	super.setValue(inValue);
    }
    
    // kliu - also need min/max
    public abstract double getMinimumValue ();
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
