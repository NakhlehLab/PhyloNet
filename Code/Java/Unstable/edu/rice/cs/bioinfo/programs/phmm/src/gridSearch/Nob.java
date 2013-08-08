package gridSearch;

import java.util.Stack;

public abstract class Nob {
    // Note: Nob is usually spelled as "Knob" but this is written in Phylonet Nob Style.

    private int g;
    private double min;
    private double max;
    protected double backupParam;

    // kliu - add in a cache for previous values
    protected Stack<Double> cache;

    public Nob(int gIn, double minIn, double maxIn) {
        this.g = gIn;
        this.min = minIn;
        this.max = maxIn;

    }

    /**
     * @return an array of sample parameter values
     * given a set number of samples and an interval
     */
    public double[] getSamples() {
        double[] samples = new double[g];
        double interval = (max - min) *1.0/ g;
        for (int iter = 0; iter < g; iter++) {
            samples[iter] = min + (interval*iter);
        }
        return samples;
    }

    /**
     * Sets a new value for the current parameter Knob
     * @param value
     */
    abstract public void set_param(double value);

    /**
     * @return the current value of the parameter knob
     */
    abstract double get_param();


    /**
     * WARNING - will throw an exception if no corresponding
     * cacheParameterValue() operation was performed previously.
     */
    abstract public void restoreParameterValue();

    // kliu - don't do a cache peek
    // unnecessary for now
}
