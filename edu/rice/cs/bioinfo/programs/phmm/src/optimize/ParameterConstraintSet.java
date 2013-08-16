/**
 * Just need some storage for constraint.
 * Immutable after creation.
 *
 * Use a BidirectionalMultimap to assign Parameter objects to a ParameterConstraintSet
 */

package optimize;

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

