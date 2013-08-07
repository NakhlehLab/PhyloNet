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
