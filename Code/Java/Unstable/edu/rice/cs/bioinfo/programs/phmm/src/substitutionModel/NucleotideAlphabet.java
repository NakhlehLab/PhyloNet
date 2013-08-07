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
