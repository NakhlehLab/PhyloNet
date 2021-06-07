/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Calculate a transition probability matrix
 * from Jukes-Cantor sequence evolution model.
 *
 * Looks good.
 *
 * Make it write once, read therafter. Since want to 
 * restrict mutability so users have to call setSubstitutionModel()
 * which calls caching functions appropriately to calculate
 * dependent cached quantities.
 */

package edu.rice.cs.bioinfo.programs.phmm.src.substitutionModel;

import edu.rice.cs.bioinfo.programs.phmm.src.util.Matrix;

public class JCSubstitutionModel implements SubstitutionModel {
    // only one parameter
    protected double substitutionRate;
    
    // only modify during initial construction
    // read only thereafter
    public JCSubstitutionModel (double inSubstitutionRate) {
	setSubstitutionRate(inSubstitutionRate);
    }

    // copy constructor
    public JCSubstitutionModel (JCSubstitutionModel jcsm) {
	setSubstitutionRate(jcsm.getSubstitutionRate());
    }

    // deep copy clone
    public SubstitutionModel deepCopyClone () {
	JCSubstitutionModel jcsm = new JCSubstitutionModel(this.substitutionRate);
	return (jcsm);
    }

    // Jukes-Cantor models have reversibility property
    public boolean checkReversible () {
	return (true);
    }

    public Alphabet getAlphabet () {
	return (NucleotideAlphabet.getClassInstance());
    }

    public double getSubstitutionRate () {
	return (substitutionRate);
    }
    
    /**
     * No more mutability after construction!
     */
    protected void setSubstitutionRate (double sr) {
	// strict!
	if (sr < 0.0) {
	    System.err.println ("ERROR: substitution rate in JCSubstitutionModel.setSubstitutionRate(double) is negative. Keeping substitution rate at 0.");
	    substitutionRate = 0.0;
	    return;
	}

	substitutionRate = sr;
    }
    
    /**
     * From Felsenstein 2001 Inferring Phylogenies
     * R[t] = e^{Qt}
     * R the substitution probability matrix
     * Q the substitution rate matrix
     * t the branch length
     */
    public double[][] calculateProbabilitiesFromRates (double time) {
	double[][] result = new double[getAlphabet().length()][getAlphabet().length()];
	for (int i = 0; i < getAlphabet().length(); i++) {
	    for (int j = 0; j < getAlphabet().length(); j++) {
		if (i == j) {
		    result[i][j] = 0.25 + 0.75 * Math.exp((-4.0 / 3.0) * time * substitutionRate);
		}
		else {
		    result[i][j] = 0.25 - 0.25 * Math.exp((-4.0 / 3.0) * time * substitutionRate);
		}
	    }
	}
	return (result);
    }

    public double[] getStationaryProbabilities () {
	double[] result = new double[getAlphabet().length()];
	for (int i = 0; i < result.length; i++) {
	    result[i] = 1.0 / ((double) getAlphabet().length());
	}
	return (result);
    }

    public String toString () {
	return ("JC substitution rate: |" + substitutionRate + "|");
    }

    // testing
    // looks fine
    public static void main (String[] args) {
	if (args.length != 2) {
	    System.err.println ("Usage: java JCSubstitutionModel <substitution rate> <branch length>");
	    System.exit(1);
	}
	JCSubstitutionModel jcsm = new JCSubstitutionModel(Double.parseDouble(args[0]));

	double[][] transitionProbabilities = jcsm.calculateProbabilitiesFromRates(Double.parseDouble(args[1]));
	
	System.out.println ("transitionProbabilities: |" + Matrix.toString(transitionProbabilities) + "|");

	double[] stationaryProbabilities = jcsm.getStationaryProbabilities();
	System.out.println ("stationaryProbabilities: |" + Matrix.toString(stationaryProbabilities) + "|");
    }

}
