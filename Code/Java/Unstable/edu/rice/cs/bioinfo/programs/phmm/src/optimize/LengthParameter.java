package optimize;

public class LengthParameter extends Parameter {
    public static final double DEFAULT_MINIMUM_BRANCH_LENGTH = MultivariateOptimizer.DEFAULT_MINIMUM_BRANCH_LENGTH;
    public static final double DEFAULT_MAXIMUM_BRANCH_LENGTH = MultivariateOptimizer.DEFAULT_MAXIMUM_BRANCH_LENGTH;

    // Use external BidirectionalMap to
    // switch between LengthParameter object and
    // edge.
    // Nice it explicitly ties between objects -
    // if multiple HiddenState objects share a Tree object,
    // we can reflect that in the map.
    //
    // Similar to t1 parameter shared across two network edges in
    // OptimizeContinuousNetworkModelParameter,
    // except that we can also share across multiple gene trees.

    public LengthParameter (String inName, double inValue) {
	super(inName, inValue);
    }

    public double getMinimumValue () {
	return (DEFAULT_MINIMUM_BRANCH_LENGTH);
    }

    public double getMaximumValue () {
	return (DEFAULT_MAXIMUM_BRANCH_LENGTH);
    }

}
