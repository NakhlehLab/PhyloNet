package gridSearch;

public abstract class Nob {
	public Nob() {
		
	}
	
	private int g;
	private double min;
	private double max;
	
	abstract public double[] getSamples();
	
	abstract public void set_param(double value);
}
