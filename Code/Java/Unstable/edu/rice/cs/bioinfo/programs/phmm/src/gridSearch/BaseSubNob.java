package gridSearch;


/**
 * This class is a Knob for the Base Subsitituion parameter
 * in the Felsenstein class.
 * @author k3kathy
 */
public class BaseSubNob extends Nob {
    private double baseRateSub = .1;	/* Default value*/

    public BaseSubNob(int gIn, double minIn, double maxIn) {
        super(gIn, minIn, maxIn);
    }

    public void set_param(double value) {
        baseRateSub = value;
    }


    public double get_param() {
        return baseRateSub;
    }

}
