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

import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.Hashtable;
import java.util.Iterator;
//import gsp.ra.Alignment;

public abstract class Alphabet {
    public static final char INDEL = '-';
    public static final char WILDCARD = 'N';

    protected Hashtable<Character,Integer> observationIndexMap;

    public abstract String getAlphabet ();

    public int length () {
	// kind of clunky
	return (getAlphabet().length());
    }

    protected String getAlphabetPlusIndel () {
	return (getAlphabet() + INDEL);
    }

    /**
     * Maintain one instance per interpreter.
     * See NucleotideAlphabet.java, for instance.
     */
    protected Alphabet () {
	initializeObservationIndexMaps();
    }

    protected void initializeObservationIndexMaps () {
	observationIndexMap = new Hashtable<Character,Integer>();

	// include indel too
	for (int i = 0; i < getAlphabetPlusIndel().length(); i++) {
	    Character C = new Character(getAlphabetPlusIndel().charAt(i));
	    Integer I = new Integer(i);
	    observationIndexMap.put(C, I);
	}
    }

    // for backwards compatibility
    public static boolean verify (String s, String alphabet) {
	return (!Pattern.matches(".*[^" + alphabet + "].*", s.toString()));
    }

    /**
     * checks if string is legal according to alphabet
     * parameterize by alphabet too since may want to modify
     */
    public boolean verify (String s, boolean includeIndelFlag, boolean includeWildcardFlag) {
	String alphabet = this.getAlphabet();
	if (includeIndelFlag) {
	    alphabet += "\\" + INDEL;
	}
	if (includeWildcardFlag) {
	    alphabet += WILDCARD;
	}
	return (verify(s, alphabet));
    }

    // public boolean verify (Alignment a, boolean includeIndelFlag, boolean includeWildcardFlag) {
    // 	Iterator<StringBuffer> sequences = a.getAllSequences();
    // 	while (sequences.hasNext()) {
    // 	    if (!verify(sequences.next().toString(), includeIndelFlag, includeWildcardFlag)) {
    // 		return (false);
    // 	    }
    // 	}
    // 	return (true);
    // }

    // public boolean verify (Alignment a) {
    // 	return (verify(a, false, false));
    // }

    /**
     * Map a character c to a unique index i in [0,|C|] for alphabet C.
     * Make this fast. Meh, not a huge deal. Later, just preprocess
     * input observation strings into strings of observation-indices.
     */
    public int getObservationSymbolIndex (char c) {
	return (observationIndexMap.get(new Character(c)).intValue());
    }


}
