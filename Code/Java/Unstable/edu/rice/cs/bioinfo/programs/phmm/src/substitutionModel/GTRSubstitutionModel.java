/**
 * Calculate a transition probability matrix
 * from sequence evolution model, e.g. GTR.
 */

package substitutionModel;

import java.util.Arrays;
import util.Matrix;
import util.Constants;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

public class GTRSubstitutionModel implements SubstitutionModel {
    protected double[][] rates;

    // need them too later
    protected double[] originalRateParameters;
    protected double[] originalFreqParameters;

    public GTRSubstitutionModel () {
	rates = new double[getAlphabet().length()][getAlphabet().length()];
    }

    public SubstitutionModel deepCopyClone () {
	GTRSubstitutionModel copy = new GTRSubstitutionModel();
	copy.originalRateParameters = Arrays.copyOf(this.originalRateParameters, this.originalRateParameters.length);
	copy.originalFreqParameters = Arrays.copyOf(this.originalFreqParameters, this.originalFreqParameters.length);

	copy.rates = new double[this.rates.length][];
	for (int i = 0; i < this.rates.length; i++) {
	    copy.rates[i] = Arrays.copyOf(this.rates[i], this.rates[i].length);
	}

	return (copy);
    }

    public double[] getStationaryProbabilities () {
	return (originalFreqParameters);
    }

    public double[] getOriginalRateParameters () {
	return (originalRateParameters);
    }

    public Alphabet getAlphabet () {
	return (NucleotideAlphabet.getClassInstance());
    }

    // GTR models have reversibility property
    public boolean checkReversible () {
	return (true);
    }

    protected int getSubstitutionParameterCount () {
	// (n choose 2) - 1
	// since all rates relative to last rate, and last rate arbitrarily set to 1.0
	return ((getAlphabet().length() * (getAlphabet().length() - 1) / 2) - 1);
    }
    
    /**
     * As it stands in NucleotideAlphabet, rateParameters is <AG> <AC> <AT> <GC> <GT>.
     * <CT> is always 1.0.
     */
    public void setSubstitutionRates (double[] rateParameters, double[] freqParameters) {
	if (rateParameters.length != getSubstitutionParameterCount()) {
	    throw (new RuntimeException("ERROR: number of rate parameters in setSubstitutionRates() should be " + getSubstitutionParameterCount() + "."));
	}

	if (freqParameters.length != getAlphabet().length()) {
	    throw (new RuntimeException("ERROR: number of stationary frequency parameters is incorrect. " + Matrix.toString(freqParameters) + " " + getAlphabet().length()));
	}

	if (!(Math.abs(Matrix.sum(freqParameters) - 1.0) <= Constants.ZERO_DELTA)) {
	    throw (new RuntimeException("ERROR: stationary frequency parameters don't sum to one. Not setting."));
	}

	originalRateParameters = rateParameters;
	originalFreqParameters = freqParameters;

	int cnt = 0;
	for (int i = 0; i < getAlphabet().length(); i++) {
	    for (int j = i + 1; j < getAlphabet().length(); j++) {
		if ((i == getAlphabet().length() - 2) && (j == getAlphabet().length() - 1)) {
		    continue;
		}
		rates[i][j] = freqParameters[j] * rateParameters[cnt];
		rates[j][i] = freqParameters[i] * rateParameters[cnt];
		cnt++;
	    }
	}

	// paranoid
	if (cnt != getSubstitutionParameterCount()) {
	    System.err.println ("ERROR: internal logic error! Count of rates not equal to calculation in setSubstitutionRates. Exiting.");
	    System.exit(1);
	}

	rates[getAlphabet().length() - 2][getAlphabet().length() - 1] = freqParameters[getAlphabet().length() - 1]; // * 1.0
	rates[getAlphabet().length() - 1][getAlphabet().length() - 2] = freqParameters[getAlphabet().length() - 2]; // * 1.0
	// sets them so that rows sum to 0
	setRateMatrixDiagonals();
	
    }

    protected void setRateMatrixDiagonals () {
	double[] rowsum = new double[rates.length];
	for (int i = 0; i < rates.length; i++) {
	    rowsum[i] = 0.0;
	    for (int j = 0; j < rates.length; j++) {
		if (i == j) {
		    continue;
		}
		rowsum[i] += rates[i][j];
	    }
	}

	for (int i = 0; i < rates.length; i++) {
	    rates[i][i] = -rowsum[i];
	}
    }
    
    /**
     * From Felsenstein 2001 Inferring Phylogenies
     * R[t] = e^{Qt}
     * R the substitution probability matrix
     * Q the substitution rate matrix
     * t the branch length
     */
    public double[][] calculateProbabilitiesFromRates (double time) {
	DoubleMatrix dm = new DoubleMatrix(rates);
	dm.mul(time);
	DoubleMatrix expdm = MatrixFunctions.expm(dm);
	return (expdm.toArray2());

	//double[][] result = new double[getAlphabet().length()][getAlphabet().length()];
    }

    public String toString () {
	return ("GTR rate matrix: |" + Matrix.toString(rates) + "|");
    }

    public static void main (String[] args) {
	GTRSubstitutionModel gsm = new GTRSubstitutionModel();
	double[] rates = new double[5];
	rates[0] = 2.0;
	rates[1] = 2.0;
	rates[2] = 2.0;
	rates[3] = 2.0;
	rates[4] = 2.0;
	double[] freqs = new double[4];
	freqs[0] = 0.25;
	freqs[1] = 0.25;
	freqs[2] = 0.25;
	freqs[3] = 0.25;
	gsm.setSubstitutionRates(rates, freqs);

	double[][] probs = gsm.calculateProbabilitiesFromRates(0.1);

	System.out.println (Matrix.toString(probs));
    }

}
