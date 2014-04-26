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

package substitutionModel;

/**
 * Eh, futz with the class hierarchy later. Geez.
 * Keep alphabet string static for efficiency reasons.
 * Only need one copy for entire program run.
 * Fix this later.
 */

import java.util.regex.Pattern;
import java.util.regex.Matcher;

public class NucleotideAlphabet extends Alphabet {
    public static final String alphabet = "AGCT";
        
    // only one instance per interpreter
    protected static final NucleotideAlphabet classInstance = new NucleotideAlphabet();

    // disallow instantiation!
    // use one instance per process pattern
    protected NucleotideAlphabet () {
	super();
    }

    public static NucleotideAlphabet getClassInstance () {
	return (classInstance);
    }

    public String getAlphabet () {
	return (alphabet);
    }

    // testing
    public static void main (String[] args) {
	if (args.length != 1) {
	    System.err.println ("Usage: java NucleotideAlphabet <string to verify>");
	    System.exit(1);
	}
	NucleotideAlphabet na = NucleotideAlphabet.getClassInstance();
	System.out.println ("alphabet: " + na.getAlphabetPlusIndel());
	System.out.println ("C " + na.getObservationSymbolIndex('C'));
	System.out.println (na.verify(args[0], true, true));
    }

}
