package runHmm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import phylogeny.EvoTree;
import phylogeny.TreeParser;
import reader.Parser;
import reader.ParserFileException;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.OpdfMap;
//import be.ac.ulg.montefiore.run.jahmm.MyHMM;
// kliu - pull in additional library support

public class runHmm {

    private static final double tolerated_error = 1e-5;		/* Sum of probabilities margin of error allowed */
	
    // stored information 
    private static String basicFileName = null;						/* Basic Info file name */
    private static String parentalTreesFileName = null;				/* Parental trees file name */
    private static String geneGenealogiesFileName = null;			/* Gene genealogies file name */
    private static String alleleSpeciesFileName = null;				/* Allele Species Mapping file name */
    

    
    // information created/built
    // kliu - index into tuples
    public static ArrayList<HiddenState> trees_states;			/* List of all states/trees */
    public static Hmm<ObservationMap> myhmm;			/* The entire HMM */
    public static double[] pi;								/* The initial pi probabilities for each state */
    public static double[][] aij;							/* The state transition probability matrix */
    private static Parser fParser;							/* The parser for all basic info and read sequences --> also calculates likelihood */
    private static int numStates = -1;						/* The number of states for the HMM */
	

    protected static void printUsage () {
	System.err.println ("Prompt-based usage: java -jar dist/lib/phmm.jar");
	System.err.println ("File-based usage: java -jar dist/lib/phmm.jar <text file with input commands>");
    }

    /**
     * kliu - add in the automatic input here.
     */
    public static void main (String[] args) throws Exception {
	if (args.length == 1) {
	    FileInputStream fis = new FileInputStream(new File(args[0]));
	    System.setIn(fis);
	    run();
	}
	else if (args.length == 0) {
	    run();
	}
	else {
	    printUsage();
	    // strict!
	    System.exit(1);
	}
    }

    /**
     * @throws Exception
     * 		- will only throw exceptions when basic Info file or Tree file is unable to be parsed correctly
     * 			since these files are essential to building an hmm and the parser
     */
    public static void run () throws Exception {
	
	int option;
	
	BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		
	option = initial(in);
	switch (option) {
	case 0: 
	    // Get the basic Info File
	    System.out.println("Input the basic file info path name:\n (note: see README for file format) \n ");
	    basicFileName = in.readLine();
			
	    // Get Trees/States Information
	    boolean getNumStates = true;
	    while (getNumStates) {
		System.out.println("Input the number of trees or states for this HMM:");
		try {
		    numStates = Integer.parseInt(in.readLine());
		    getNumStates = false;
		} catch (NumberFormatException e) {
		    System.out.println( e + "\nError: Not a number! Try again!");
		}
	    }
			
	    System.out.println("\nInput the parental trees file path name:\n (note: see README for file format) \n");
	    parentalTreesFileName = in.readLine();
			
	    System.out.println("\nInput the gene genealogies file path name:\n (note: see README for file format) \n");
	    geneGenealogiesFileName = in.readLine();
	    
	    System.out.println("\nInput the Alleles to Species Mapping file path name: \n (note: see README for file format \n");
	    alleleSpeciesFileName = in.readLine();

	    // Get the Pi probabilities array
	    getPiInfo(in);
	    
	    // Get transition Matrix
	    getAij(in);
			
	    // Now build the trees and get the tree mapping to integers
	    buildTrees();
	    
	    //Reading in Basic Info file and store information
	    buildParser();
			
	    // Building Allele to Species map
	    buildAlleleSpeciesMap();
	    
	    //Build HMM
	    buildHMM();
		
	    // strict!
	    if (!verifyInputs()) {
		System.err.println ("ERROR: inputs are invalid. Please correct and try again.");
		System.exit(1);
	    }
	    
	    break;
	    
	case 1: 
	    System.out.println("Get existing model");
	    break;
		
	case 2: 
	    System.exit(0);
			
	default:
	    System.exit(-1);
	}
		
		
		
	//Operate mode
	boolean keepOperate = true;
	while (keepOperate) {
	    option = operate(in);
	    switch(option) {
	    case 0:
		    // RUN VITERBI
			System.out.print("\n");
			// get name of output file
			System.out.println("Path to your output file: ");
			String outputfile = in.readLine();
					
			// kliu - only Viterbi algorithm appears to be implemented using Jahmm
			// don't seem to have forward/backward implementation yet
			ArrayList<ObservationMap> obsSeq = getObs(in);
			
			myhmm.saveMostLikelyStateSequence(obsSeq, outputfile);
					
			break;
	    case 1:
		    // RUN BAUM WELCH
			
			// Testing purposes
//			ArrayList<ObservationMap> oSeq = listOfOseq.get(0);
//			System.out.println("\n ---------------------\nThe HMM Emission Probabilities");
//			for (int i = 0; i < trees_states.size(); i++) {
//				System.out.println("State : " + i);
//				System.out.println("Opdf : {");
//				for (int j = 0; j < oSeq.size(); j++) {
//					System.out.print(oSeq.get(j) + ":" + myhmm.getOpdf(i).probability(oSeq.get(j)) + " , ");
//				}
//				System.out.print("}\n\n");
//			}
			
	    	// END TEST PURPOSES
				
				
			BaumWelchScaledLearner bwl = new BaumWelchScaledLearner();
			
			Hmm<ObservationMap> learnedhmm = bwl.learn(myhmm,getObs(in), 9);
			myhmm = learnedhmm;
			learnedhmm = null; // for garbage collection
			System.out.println("\n----------------------------\nThe LEARNED HMM IS: ");
			System.out.println(myhmm + "\n\n");
			
			break;
		
		
	    case 2:
	    // EXIT
			keepOperate = false;
			System.exit(0);
			break;
	    default: 
			keepOperate = false;
			System.exit(-1);
	    }
	}
    }
    
    /**
     * Additional verification step.
     */
    static protected boolean verifyInputs () {
	if (!verifyGeneGenealogyTaxa()) {
	    System.err.println ("ERROR: unable to verify gene genealogy inputs.");
	    return (false);
	}
	
	// check allele-to-species mapping
	if (!verifyParentalTreeTaxa()) {
	    System.err.println ("ERROR: unable to verify parental tree inputs.");
	    return (false);
	}
	
	// others as needed

	return (true);
    }
    
    /**
     * Warning - constructs a new array each time.
     */
    static protected String[] getSortedGeneGenealogyTaxa () {
	// array copy - don't lose original list order
	String[] sortedGeneGenealogyTaxa = new String[fParser.getGeneGenealogyTaxa().size()];
	sortedGeneGenealogyTaxa = fParser.getGeneGenealogyTaxa().toArray(sortedGeneGenealogyTaxa);
	Arrays.sort(sortedGeneGenealogyTaxa);
	return (sortedGeneGenealogyTaxa);
    }

    /**
     * Use allele-to-species mapping.
     * Set keyFlag to true to get keys (allele taxa), or false to get values (species taxa).
     * Warning - constructs a new array each time.
     */
    static protected String[] getSortedTaxaFromAlleleSpeciesMapping (boolean keyFlag) {
	HashMap<String,String> alleleSpeciesMap = fParser.getAlleleSpeciesMap();
	// no official list of species names - meh
        Set<String> set = keyFlag ? new HashSet<String>(alleleSpeciesMap.keySet()) : new HashSet<String>(alleleSpeciesMap.values());
	String[] sortedParentalTreeTaxa = new String[set.size()];
	sortedParentalTreeTaxa = set.toArray(sortedParentalTreeTaxa);
	Arrays.sort(sortedParentalTreeTaxa);
	return (sortedParentalTreeTaxa);
    }

    static protected boolean verifyGeneGenealogyTaxa () {
	// list of taxa in allele-species-mapping and in basic-info-file must be consistent with 
	// gene genealogies
	if (!verifyTaxa(getSortedGeneGenealogyTaxa(), true)) {
	    System.err.println ("ERROR: list of taxa in basic-info-file is not consistent with gene genealogies.");
	    return (false);
	}

	if (!verifyTaxa(getSortedTaxaFromAlleleSpeciesMapping(true), true)) {
	    System.err.println ("ERROR: list of taxa in allele-species-mapping is not consistent with gene genealogies.");
	    return (false);
	}

	return (true);
    }

    static protected boolean verifyParentalTreeTaxa () {
	if (!verifyTaxa(getSortedTaxaFromAlleleSpeciesMapping(false), false)) {
	    System.err.println ("ERROR: list of taxa in basic-info-file is not consistent with parental trees.");
	    return (false);
	}

	return (true);
    }

    /** 
     * Set geneGenealogyTaxaFlag to true to check gene genealogy taxa, set to false
     * to check parental tree taxa.
     */
    static protected boolean verifyTaxa (String[] sortedReferenceTaxa, boolean geneGenealogyTaxaFlag) {
	for (HiddenState hiddenState : trees_states) {
	    String[] sortedCheckTaxa = geneGenealogyTaxaFlag ? 
		hiddenState.getGeneGenealogy().getTaxa() : 
		hiddenState.getParentalTree().getTaxa();
	    Arrays.sort(sortedCheckTaxa);
	    if (!Arrays.equals(sortedReferenceTaxa, sortedCheckTaxa)) {
		return (false);
	    }
	}

	return (true);
    }


    /**
     * Initial step and options:
     * a. build new model
     * b. load a pre-existing model
     * c. exit
     * 
     * @return the option
     */
    static private int initial(BufferedReader in) {
	boolean init = true;
	int option = -1;
	while (init) {
	    System.out.println("Initial mode:");
	    System.out.println("0) Build a new model.");
	    System.out.println("1) Load a pre-existing model.");
	    System.out.println("2) Exit.");
	    System.out.println("Choose an option: ");
	    option = getOption(3, in);
	    if (option != -1) init = false;
	}
	return option;
    }
	
	
	
    /**
     * Operate on Model:
     * a. Read in sequence files
     * b. Use Viterbi's on existing model
     * c. exit
     * 
     * @return the option
     */
    private static int operate(BufferedReader in) {
	boolean operate = true;
	int option = -1;
	while (operate) {
	    System.out.println("Operate mode:");
	    System.out.println("0) Process an observation sequence");
	    System.out.println("1) Learn Model using Baum Welch.");
	    System.out.println("2) Exit.");
	    System.out.println("Choose an option: ");
	    option = getOption(3, in);
	    if (option != -1) operate = false;
	}
	return option;
    }
	
	

	

    /**
     * Read and Get Observation
     */
    private static ArrayList<ObservationMap> getObs(BufferedReader in) {
	try {
	    System.out.println("Input observation sequence file path name : ");
	    String filename = in.readLine();
	    fParser.parseMe(filename, myhmm);
	    return fParser.getObs();
	} catch (Exception e) {
	    e.printStackTrace();
	}
	return null;
		
    }
	
	
    /**
     * Helper function
     * Given the number of options, this method will read in the user
     * choice while catching illegal inputs and output the user choice option
     * 
     * @param numOptions int - the number options available
     * @param option int - the option the user inputted
     * @return int the option number
     */
    private static int getOption(int numOptions, BufferedReader in) {
	try {
	    int option = Integer.parseInt(in.readLine());
	    if ((option >= 0) && (option < numOptions)) {
		return option;
	    }
	    else {
		System.out.println("Your answer is not an option!");
		return -1;
	    }
	} catch (NumberFormatException e) {
	    System.out.print("Error: not a number!");
	    return -1;
	} catch (IOException e) {
	    e.printStackTrace();
	    return -1;
	}
    }
	
    // kliu - no reason to put it in a separate class
    public static Hmm<ObservationMap> buildMyHmm(List<HiddenState> hiddenStates, double[] pi, double[][] a) {		
	if ((pi.length != a.length) || (hiddenStates.size() != pi.length)) {
	    throw new IllegalArgumentException("The dimension of the initial probability array, the transition probability matrix number of rows, and the number of hidden states are unequal.");
	}
		
	int nbStates = pi.length;

	// kliu testing
	//System.out.println ("foo: |" + nbObs + "|");
	
	ArrayList<OpdfMap> opdfs = new ArrayList<OpdfMap>(nbStates);	
	//ArrayList<Opdf<ObservationInteger>> opdfs = new ArrayList<Opdf<ObservationInteger>>(nbStates);
		
	for (int i = 0; i < nbStates; i++) {
	    opdfs.add(new OpdfMap(hiddenStates.get(i)));
	}
		
	// kliu - needs to not be cached - just need to map back from Opdf
	// to the hidden state that's associated with it -> 
	Hmm<ObservationMap> newHmm = new Hmm<ObservationMap>(pi, a, opdfs);
		
	return newHmm;		
    }

	
    /**
     * Build HMM
     */
    private static void buildHMM() {
	System.out.println("\n\nNow building HMM . . .");
	myhmm = buildMyHmm(trees_states, pi, aij);
	//myhmm = MyHMM.buildMyHmm(trees_states, pi, aij);
	System.out.println("\nhey my hmm: \n " + myhmm);
    }

//(int) Math.pow(fParser.getAlphabet().size(),fParser.getNumSeq()
	
    /**
     * Read and parse basic info file
     * to build the Parser for reading sequences
     * @throws Exception 
     */
    private static void buildParser() throws Exception {
		if(basicFileName == null) {
			throw new ParserFileException("Cannot read Basic Info File!");
		}
	    System.out.println("\nNow reading and saving Basic Info for parser . . .");
	    fParser = new Parser(basicFileName);
	    fParser.setTrees(trees_states);	
    }
	
    
    /**
     * Reads in and parse the alleles to species map
     * @throws Exception
     */
    private static void buildAlleleSpeciesMap() throws Exception {
    	if (alleleSpeciesFileName == null) {
    	    throw new ParserFileException("Cannot read Trees file!");
    	}
    	
    	fParser.parseAlleleSpecies (alleleSpeciesFileName);
	
	// kliu - paranoid
	if ((fParser.getAlleleSpeciesMap() == null) ||
	    fParser.getAlleleSpeciesMap().isEmpty()) {
	    // strict!
	    System.err.println ("ERROR: allele-to-species map empty or incorrectly formatted. Aborting.");
	    System.exit(1);
	}

	// now set references in all hidden states
	for (HiddenState hiddenState : trees_states) {
	    hiddenState.setAlleleToSpeciesMapping(fParser.getAlleleSpeciesMap());
	}

    	//Testing Purposes//
//    	HashMap<String,String> amap = fParser.getAlleleSpeciesMap();
//    	Set<String> keyset = amap.keySet();
//    	for (String i : keyset) {
//    		System.out.println(i + " : " + amap.get(i));
//    	}
    	// Testing purposes//
    }
	
    /**
     * Read and parse trees file
     * @throws Exception 
     */
    private static void buildTrees() throws Exception {
	if (parentalTreesFileName == null) {
	    throw new ParserFileException("Cannot read Trees file!");
	}

	trees_states = new ArrayList<HiddenState>();

	System.out.println("\nNow building trees . . .");
	BufferedReader ptreesbr = new BufferedReader(new FileReader(parentalTreesFileName));
	TreeParser ptp = new TreeParser(ptreesbr);
	ArrayList<EvoTree> parentalTrees = ptp.nexusFileTreeNames(parentalTreesFileName);
	ptreesbr.close();

	// kliu - indexing is by (parentalTree, geneGenealogy) appearance order according to the following:
	for (EvoTree parentalTree : parentalTrees) {
	    // kliu - cheap hack to clone all gene genealogies across
	    // each parental tree
	    //
	    // for now, don't share gene genealogies across parental trees (e.g. branch lengths)
	    BufferedReader ggbr = new BufferedReader(new FileReader(geneGenealogiesFileName));
	    TreeParser gtp = new TreeParser(ggbr);
	    ArrayList<EvoTree> geneGenealogies = gtp.nexusFileTreeNames(geneGenealogiesFileName);
	    ggbr.close();

	    for (EvoTree geneGenealogy : geneGenealogies) {
		// kliu - meh - parse allele-to-species mapping later and add in references here
		trees_states.add(new HiddenState(parentalTree, geneGenealogy, null));
	    }
	}
	
	// ------- >Testing purposes
	//	        System.out.println("Trees read in and built:");
	//	        for (int i = 0; i < trees_states.size(); i++) {
	//	        	System.out.println((trees_states.get(i)));
	//	        }
	// <------Testing purposes

	if (numStates != trees_states.size()) {
	    throw new ParserFileException("Error: the number of trees read is not the same as the number of states previously inputted!");
	}	
    }
	
	
	
	
    /**
     * Get State Transition Matrix Aij
     * @throws IOException - rare
     */
    private static void getAij(BufferedReader in) throws IOException {
	boolean getA = true;
	boolean getChoice = true;
		
	while (getA) {			
	    // Getting transition matrix A
	    System.out.println("\nThe Transition of States/Trees Matrix Aij. \n(Note: Make sure the order in the matrix correspond to the same order of trees in the given tree input file.");
	    while (getChoice) {
		System.out.println("Choose an option: \n(a) Read transition matrix by file.\n(b) Input transition matrix manually.");
		String choice = in.readLine().toLowerCase();
		if (choice.equals("a")) {
		    try {
			System.out.println("Please see README for .matrix file format.");
			System.out.println("Input .matrix file path name:");
			String filename = in.readLine();
			BufferedReader filebr = new BufferedReader(new FileReader(filename));
			aij = readAij(numStates, filebr);
			filebr.close();
			getChoice = false;
		    } catch (Exception e) {
			System.out.println(e);
			System.out.println("Please try again. \n");
		    }
		}
		else if (choice.equals("b")) {
		    try {
			System.out.println("Example transition for 3 trees (inputted 3 times): \n   .3 .2 .5 \n   .2 .4 .4 \n   .5 .4 .1 \n");
			System.out.println("Input transition Aij matrix:");
			aij = readAij(numStates, in);
			getChoice = false;
		    } catch (Exception e) {
			System.out.println(e);
			System.out.println("Please try again.");
		    }
		}
	    }
	    getA = false;
	}
    }
	
	
	
    /**
     * Get Initial Pi Probabilities array
     * @throws IOException - rare
     */
    private static void getPiInfo(BufferedReader in) throws IOException {
	boolean getPi = true;
		
	while (getPi) {
			
	    // Getting pi
	    System.out.println("\nThe Pi or Initial Probabilities. \n(Note: Make sure pi probabilities correspond to the same order of trees in the trees input file.) \n Sample input for 4 trees/states: .2 .3 .1 .4");
	    System.out.println("Input Initial Pi probabilities:");
	    String[] piString;
	    piString = in.readLine().split(" ");

				
	    // Error Check: check to see if number of states and number of pi probabilities are the same
	    if (numStates == piString.length) {
		getPi = false;
					
		// make pi probabilities array
		pi = new double[numStates];
		double sum = 0;
		for (int i=0; i < piString.length; i++) {
		    try {
			pi[i] = Double.parseDouble(piString[i]);
			sum += pi[i];
		    } 
		    catch (NumberFormatException e) { 
			getPi = true;
			System.out.println("Number format error: Cannot convert inputed pi number : " + piString[i] + " to a double."); 
		    }
		}
		if (Math.abs(1.0 - sum) > tolerated_error) {
		    System.out.println("Error: the sum of the pi probabilities did not sum up to 1.0");
		    getPi = true;
		}
					
	    } else System.out.println("Error: number of pi probability inputs did not match number of states! ");

	}
    }
	
	
	
    /**
     * Read in Transition matrix
     * @param numStates - the number of states/trees
     * @param br - a buffered reader
     * @returns a double[][] transition of states/trees matrix probabilities
     * @throws Exception - Either an IO exception if an IO error occurred or a ParserFileException that signals an incorrect matrix input.
     */
    public static double[][] readAij(int numStates, BufferedReader br) throws Exception {
	String[] row;
	double[] newRow;
	double[][] aij = new double[numStates][numStates];
		
	for (int i = 0; i < numStates; i++) {
	    newRow = new double[numStates];
	    row = br.readLine().split(" ");
	    double sum = 0;
	    if (row.length != numStates) throw new ParserFileException("Error: Number of entries in matrix is incorrect!");
	    for (int j = 0; j < row.length; j++) {
		newRow[j] = Double.parseDouble(row[j]);
		sum += newRow[j];
	    }
			
	    // Check to make sure the sum of the entire row sums up to 1.0
	    if (Math.abs(1.0 - sum) > tolerated_error) {
		throw new ParserFileException("Error: the sum of row " + i + " did not sum up to 1.0.");
	    }
			
	    aij[i] = newRow;
	}
		
	// Check to make sure the sum of each column sums up to 1.0
	for (int j = 0; j < numStates; j++) {
	    double sum = 0;
	    for (int i = 0; i < numStates; i++) {
		sum+= aij[i][j];
	    }
	    if (Math.abs(1.0 - sum) > tolerated_error) {
		throw new ParserFileException("Error: the sum of col " + j + " did not sum up to 1.0.");
	    }
	}
	return aij;
    }
	
}
