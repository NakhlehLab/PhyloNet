package gridSearch;

public abstract class Nob {
	public Nob() {
		
	}
	
	private int g;
	private double min;
	private double max;
	
	/*
	 * returns an array of sample parameter values
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
	
	abstract public void set_param(double value);
}
