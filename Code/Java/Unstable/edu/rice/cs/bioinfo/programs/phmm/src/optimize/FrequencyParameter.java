package optimize;

public class FrequencyParameter extends Parameter {
    public static final double DEFAULT_MINIMUM_PROBABILITY = MultivariateOptimizer.DEFAULT_MINIMUM_PROBABILITY;
    public static final double DEFAULT_MAXIMUM_PROBABILITY = MultivariateOptimizer.DEFAULT_MAXIMUM_PROBABILITY;

    public FrequencyParameter (String inName, double inValue) {
    super(inName, inValue);
    }

    public double getMinimumValue () {
    return (DEFAULT_MINIMUM_PROBABILITY);
    }

    public double getMaximumValue () {
    return (DEFAULT_MAXIMUM_PROBABILITY);
    }

}
