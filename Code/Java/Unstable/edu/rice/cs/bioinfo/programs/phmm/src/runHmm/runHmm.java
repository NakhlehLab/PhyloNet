package runHmm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.Vector;

import optimize.MultivariateOptimizer;
import gridSearch.GridSearchAlgorithm;
import phylogeny.EvoTree;
import phylogeny.TreeParser;
import reader.Parser;
import reader.ParserFileException;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.OpdfMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.library.programming.BijectiveHashtable;

public class runHmm {

    protected static final double tolerated_error = 1e-5;		/* Sum of probabilities margin of error allowed */
	
    // stored information 
    protected String basicFileName = null;						/* Basic Info file name */
    protected String parentalTreesFileName = null;				/* Parental trees file name */
    protected String geneGenealogiesFileName = null;			/* Gene genealogies file name */
    protected String alleleSpeciesFileName = null;				/* Allele Species Mapping file name */
    

    
    // information created/built
    // kliu - index into tuples
    protected ArrayList<HiddenState> trees_states;			/* List of all states/trees */
    protected Hmm<ObservationMap> myhmm;			/* The entire HMM */

    // argh - in lieu of worrying about EvoTree.equals() method
    // Maintain equivalence classes among hidden states based on shared parental tree.
    protected Map<Network<Double>,Set<HiddenState>> parentalTreeClasses;
    // Since it's bijective, make it reversible
    protected BijectiveHashtable<String,Network<Double>> parentalTreeNameMap;

    protected TransitionProbabilityParameters transitionProbabilityParameters;

    // kliu - neither are necessary anymore
    //    public static double[] pi;								/* The initial pi probabilities for each state */    
    //    public static double[][] aij;							/* The state transition probability matrix */
    protected Parser fParser;							/* The parser for all basic info and read sequences --> also calculates likelihood */
    protected int numStates = -1;						/* The number of states for the HMM */
	
    protected static void printUsage () {
	System.err.println ("Prompt-based usage: java -jar dist/lib/phmm.jar");
	System.err.println ("File-based usage: java -jar dist/lib/phmm.jar <text file with input commands>");
    }

    /**
     * kliu - add in the automatic input here.
     */
    public static void main (String[] args) throws Exception {
	runHmm rh = new runHmm();
	if (args.length == 1) {
	    FileInputStream fis = new FileInputStream(new File(args[0]));
	    System.setIn(fis);
	    rh.run();
	}
	else if (args.length == 0) {
	    rh.run();
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
    public void run () throws Exception {
	
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

	    readTransitionProbabilityParameters(in);

	    // Get the Pi probabilities array
	    //getPiInfo(in);
	    
	    // Get transition probability parameters.
	    // See writeup for details.
	    //getAij(in);
			
	    // Now build the trees and get the tree mapping to integers
	    buildTrees();
	    
	    //Reading in Basic Info file and store information
	    buildParser();
			
	    // Building Allele to Species map
	    buildAlleleSpeciesMap();
	    
	    //Build HMM
	    buildInitialHMM();
		
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

			// Allow users to read in New Sequence Observation or re-use a
            // previously read in Sequence
			ArrayList<ObservationMap> obsSeq;
			int seqChoice = sequenceChoice(in);
            if (seqChoice == 1) {
                obsSeq = getObs(in);
            } else {
                obsSeq = fParser.getObs();
            }


			double viterbiLLH = myhmm.saveMostLikelyStateSequence(obsSeq, outputfile);
			
			System.out.println ("Input HMM Viterbi log likelihood: |" + viterbiLLH + "|");

			System.out.println ("Computing input HMM log likelihood for input sequences... ");
			double llh = myhmm.lnProbability(obsSeq);
			System.out.println ("Computing input HMM log likelihood for input sequences DONE.");
			
			System.out.println ("Input HMM log likelihood: |" + llh + "|");
	
			break;
	    case 1:
		    // RUN BAUM WELCH
	        // Modifies the current Hmm;

	        System.out.println("\n");

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

	        // Allow users to read in New Sequence Observation or re-use a
	        // previously read in Sequence
	        int seqChoice2 = sequenceChoice(in);
	        ArrayList<ObservationMap> obsSequence;
	        if (seqChoice2 == 1) {
	            obsSequence = getObs(in);
	        } else {
	            obsSequence = fParser.getObs();
	        }

			BaumWelchScaledLearner bwl = new BaumWelchScaledLearner();

			Hmm<ObservationMap> learnedhmm = bwl.learn(myhmm, obsSequence, 9);
			myhmm = learnedhmm;
			learnedhmm = null; // for garbage collection
			System.out.println("\n----------------------------\nThe LEARNED HMM IS: ");
			System.out.println(myhmm + "\n\n");
			
			break;
		
		
	    case 2:
	        // RUN GRID SEARCH
	        // Modifies the current hmm;

	        System.out.println("\n");
	        // Allow users to read in New Sequence Observation or re-use a
            // previously read in Sequence
            int seqChoice3 = sequenceChoice(in);
            ArrayList<ObservationMap> seqObs;
            if (seqChoice3 == 1) {
                seqObs = getObs(in);
            } else {
                seqObs = fParser.getObs();
            }

            //Get the max and mins and g samples of each parameter
            try {
                System.out.println("Input the number of samples for tree branch lengths: ");
                int gBranchIn = Integer.parseInt(in.readLine());

                System.out.println("Input the number of samples for Recombination Frequency: ");
                int gRecombinationIn = Integer.parseInt(in.readLine());

                System.out.println("Input the number of samples for Hybridization Frequency: ");
                int gHybridizationIn = Integer.parseInt(in.readLine());

                System.out.println("Input the number of samples for Felsenstein Base Substitution Rate: ");
                int gBaseSubIn = Integer.parseInt(in.readLine());

                System.out.println("Input the minimum double value for all tree branch lengths: ");
                double branchMinIn = Double.parseDouble(in.readLine());

                System.out.println("Input the maximum double value for all tree branch lengths: ");
                double branchMaxIn = Double.parseDouble(in.readLine());

                System.out.println("Input the minimum double value for Recombination frequency: ");
                double recombinationMinIn = Double.parseDouble(in.readLine());

                System.out.println("Input the maximum double value for Recombination frequency: ");
                double recombinationMaxIn = Double.parseDouble(in.readLine());

                System.out.println("Input the minimum double value for Hybridization frequency: ");
                double hybridizationMinIn = Double.parseDouble(in.readLine());

                System.out.println("Input the maximum double value for Hybridization frequency: ");
                double hybridizationMaxIn = Double.parseDouble(in.readLine());

                System.out.println("Input the minimum double value for Felsenstein base substitution rate: ");
                double baseSubMinIn = Double.parseDouble(in.readLine());

                System.out.println("Input the maximum double value for Felsenstein base substitution rate: ");
                double baseSubMaxIn = Double.parseDouble(in.readLine());

                GridSearchAlgorithm gsa =
                        new GridSearchAlgorithm(gBranchIn, gRecombinationIn,
                        gHybridizationIn, gBaseSubIn, branchMinIn, branchMaxIn,
                        recombinationMinIn, recombinationMaxIn, hybridizationMinIn,
                        hybridizationMaxIn, baseSubMinIn, baseSubMaxIn);

                gsa.runGridSearch(seqObs, myhmm, transitionProbabilityParameters,
                        trees_states, parentalTreeClasses);

            }
            catch (Exception e) {
                System.err.println("Input error : " + e);
                break;
            }

	        break;
	    case 3:
		System.out.println("\n");
		runMultivariateOptimizer(in);
		break;
	    case 4:
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

    protected void runMultivariateOptimizer (BufferedReader in) {
	try {
	    int sc = sequenceChoice(in);
	    ArrayList<ObservationMap> obsSequence;
	    if (sc == 1) {
		obsSequence = getObs(in);
	    } 
	    else {
		obsSequence = fParser.getObs();
	    }

	    System.out.println("Input length-parameter-to-edge map filename: ");
	    String inputLengthParameterToEdgeMapFilename = in.readLine();
	    System.out.println("Input length-parameter-set-constraints filename: ");
	    String inputLengthParameterSetConstraintsFilename = in.readLine();
	    System.out.println("Output Viterbi-optimal hidden state sequence file: ");
	    String viterbiHiddenStateSequenceFilename = in.readLine();
	    System.out.println("Output model likelihoods file: ");
	    String modelLikelihoodsFilename = in.readLine();

	    MultivariateOptimizer multivariateOptimizer = new MultivariateOptimizer(myhmm,
										    this,
										    trees_states,
										    transitionProbabilityParameters,
										    parentalTreeClasses,
										    parentalTreeNameMap,
										    obsSequence,
										    inputLengthParameterToEdgeMapFilename,
										    inputLengthParameterSetConstraintsFilename
										    );

	    System.out.println ("Optimizing PhyloNet-HMM parameters... ");
	    multivariateOptimizer.optimize();
	    System.out.println ("Optimizing PhyloNet-HMM parameters DONE. ");

	    System.out.println ("Saving Viterbi-optimal hidden state sequence... "); 
	    double viterbiLLH = myhmm.saveMostLikelyStateSequence(obsSequence, viterbiHiddenStateSequenceFilename);
	    System.out.println ("Saving Viterbi-optimal hidden state sequence DONE."); 		

	    System.out.println ("Computing final model log likelihood... ");
	    double finalLLH = myhmm.lnProbability(obsSequence);
	    System.out.println ("Computing final model log likelihood DONE. ");

	    System.out.println ("Optimized model's final log likelihood: |" + finalLLH + "|");
	    System.out.println ("Optimized model's Viterbi-optimal log likelihood: |" + viterbiLLH + "|");

	    System.out.println ("Saving log likelihoods... ");
	    BufferedWriter bw = new BufferedWriter(new FileWriter(modelLikelihoodsFilename));
	    bw.write(Double.toString(finalLLH)); bw.newLine();
	    bw.write(Double.toString(viterbiLLH)); bw.newLine();
	    bw.flush();
	    bw.close();
	    System.out.println ("Saving log likelihoods DONE. ");
	}
	catch (IOException ioe) {
	    System.err.println (ioe);
	    ioe.printStackTrace();
	}
    }
    
    /**
     * Additional verification step.
     */
    protected boolean verifyInputs () {
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
    protected String[] getSortedGeneGenealogyTaxa () {
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
    protected String[] getSortedTaxaFromAlleleSpeciesMapping (boolean keyFlag) {
	HashMap<String,String> alleleSpeciesMap = fParser.getAlleleSpeciesMap();
	// no official list of species names - meh
        Set<String> set = keyFlag ? new HashSet<String>(alleleSpeciesMap.keySet()) : new HashSet<String>(alleleSpeciesMap.values());
	String[] sortedParentalTreeTaxa = new String[set.size()];
	sortedParentalTreeTaxa = set.toArray(sortedParentalTreeTaxa);
	Arrays.sort(sortedParentalTreeTaxa);
	return (sortedParentalTreeTaxa);
    }

    protected boolean verifyGeneGenealogyTaxa () {
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

    protected boolean verifyParentalTreeTaxa () {
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
    protected boolean verifyTaxa (String[] sortedReferenceTaxa, boolean geneGenealogyTaxaFlag) {
	for (HiddenState hiddenState : trees_states) {
	    String[] sortedCheckTaxa = geneGenealogyTaxaFlag ? 
		hiddenState.getGeneGenealogy().getLeaves() : 
		getTaxa(hiddenState.getParentalTree());
	    Arrays.sort(sortedCheckTaxa);
	    if (!Arrays.equals(sortedReferenceTaxa, sortedCheckTaxa)) {
		return (false);
	    }
	}

	return (true);
    }

    /**
     * Check for duplicate taxa in PhyloNet tree object.
     */
    protected boolean hasDuplicateNames (Tree tree) {
	String[] names = tree.getLeaves();
	HashSet<String> hs = new HashSet<String>();
	for (String name : names) {
	    if (hs.contains(name)) {
		return (true);
	    }
	    hs.add(name);
	}
	return (false);
    }

    /**
     * Get names of taxa from a PhyloNet network object.
     */
    protected String[] getTaxa (Network<Double> network) {
	Vector<String> names = new Vector<String>();
	for (NetNode<Double> leaf : network.getLeaves()) {
	    names.add(leaf.getName());
	}
	String[] result = new String[names.size()];
	return (names.toArray(result));
    }

    /**
     * Initial step and options:
     * a. build new model
     * b. load a pre-existing model
     * c. exit
     * 
     * @return the option
     */
    protected int initial(BufferedReader in) {
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
    protected int operate(BufferedReader in) {
	boolean operate = true;
	int option = -1;
	while (operate) {
	    System.out.println("Operate mode:");
	    System.out.println("0) Run Viterbi");
	    System.out.println("1) Learn Model using Baum Welch.");
	    System.out.println("2) Learn Model using Grid Search");
	    System.out.println("3) Learn Model using a multivariate optimization heuristic that incorporates Brent's method");
	    System.out.println("4) Exit");
	    System.out.println("Choose an option: ");
	    option = getOption(5, in);
	    if (option != -1) operate = false;
	}
	return option;
    }
	
	


    /**
     * Re-Use sequence or Insert a new Sequence?
     *
     * @return the option
     */
    protected int sequenceChoice(BufferedReader in) {
    boolean init = true;
    int option = -1;
    while (init) {
        System.out.println("Which observation sequence would you like to use?");
        System.out.println("0) Reuse previously read sequence.");
        System.out.println("1) Load a new observation sequence.");
        System.out.println("Choose an option: ");
        option = getOption(2, in);
        if (option != -1) init = false;
    }
    return option;
    }


    /**
     * Read and Get Observation
     */
    protected ArrayList<ObservationMap> getObs(BufferedReader in) {
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
    protected int getOption(int numOptions, BufferedReader in) {
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
    public Hmm<ObservationMap> buildMyHmm(List<HiddenState> hiddenStates, double[] pi, double[][] a) {		
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
     * Initial state probabilities.
     * \pi(s_i) = z(s_i) / (\sum_{s \in S} z(s)
     * where z(s) is P[g(s_j)|T(s_j), c_{T(s_k)}]
     *
     * See writeup for details.
     * 
     * Also call this after changing parental tree branch lengths.
     */
    protected double[] calculatePi () {
	double[] pi = new double[trees_states.size()];
	double norm = 0.0;
	for (int i = 0; i < trees_states.size(); i++) {
	    HiddenState hiddenState = trees_states.get(i);
	    pi[i] = hiddenState.calculateProbabilityOfGeneGenealogyInParentalTree();
	    norm += pi[i];
	}
	
	// barf if norm is basically zero
	if (norm < tolerated_error) {
	    System.err.println ("ERROR: sum of initial hidden state distribution is zero in calculateInitialPi(). Returning null to signal error.");
	    return (null);
	}

	// normalize
	for (int i = 0; i < pi.length; i++) {
	    pi[i] /= norm;
	}

	return (pi);
    }

    protected boolean checkSameParentalClass (HiddenState si, HiddenState sj) {
	Set<HiddenState> sic = parentalTreeClasses.get(si.getParentalTree());
	return (sic.contains(sj));
    }

    /**
     * Calculate initial transition probability matrix a_{ij}.
     * See revised writeup for details.
     *
     * Call this after changing parental tree branch lengths
     */
    protected double[][] calculateAij () {
	double[][] a = new double[trees_states.size()][trees_states.size()];
	for (int i = 0; i < a.length; i++) {
	    HiddenState si = trees_states.get(i);
	    double totalNonSelfTransitionProbabilities = 0.0;
	    for (int j = 0 ; j < a[i].length; j++) {
		// set self-transition probability at the end
		if (i == j) {
		    continue;
		}
		
		HiddenState sj = trees_states.get(j);
		a[i][j] = sj.calculateProbabilityOfGeneGenealogyInParentalTree();
		if (checkSameParentalClass(si, sj)) {
		    a[i][j] *= transitionProbabilityParameters.getRecombinationFrequency();
		}
		else {
		    a[i][j] *= transitionProbabilityParameters.getHybridizationFrequency();
		}

		totalNonSelfTransitionProbabilities += a[i][j];
	    }
	    
	    // now set self-transition probability
	    a[i][i] = 1.0 - totalNonSelfTransitionProbabilities;
	}
	
	// strict!
	if (!verifyAij(a)) {
	    System.err.println ("ERROR: verifyAij() failed. Returning null to signal error.");
	    return (null);
	}

	return (a);
    }

    /**
     * Given the current model parameter values, 
     * calculate the maximum possible value for a frequency parameter.
     * Anything greater will cause a self-transition probability to go negative.
     *
     * Hmm... need to expose this externally.
     */
    public double calculateMaximumFrequencyParameter (TransitionProbabilityParameters.ParameterChoice parameterChoice) {
	double minmax = -1.0;
	boolean initializedFlag = false; // bleh

	// strict!
	if (trees_states.size() <= 1) {
	    throw (new RuntimeException("ERROR: not enough hidden states in runHmm.calculateMaximumFrequencyParameter(...)."));
	}

	// take min max value over all states
	for (int i = 0; i < trees_states.size(); i++) {
	    HiddenState si = trees_states.get(i);
	    double sumCoalescentContributionRecombination = 0.0;
	    double sumCoalescentContributionHybridization = 0.0;

	    // compute max for a state si
	    for (int j = 0 ; j < trees_states.size(); j++) {
		// set self-transition probability at the end
		if (i == j) {
		    continue;
		}
		
		HiddenState sj = trees_states.get(j);
		if (checkSameParentalClass(si, sj)) {
		    sumCoalescentContributionRecombination += sj.calculateProbabilityOfGeneGenealogyInParentalTree();
		}
		else {
		    sumCoalescentContributionHybridization += sj.calculateProbabilityOfGeneGenealogyInParentalTree();
		}
	    }

	    double max;
	    switch (parameterChoice) {
	    case RECOMBINATION_FREQUENCY:
		max = (1 - sumCoalescentContributionHybridization * transitionProbabilityParameters.get(TransitionProbabilityParameters.ParameterChoice.HYBRIDIZATION_FREQUENCY)) / sumCoalescentContributionRecombination;
		break;
	    case HYBRIDIZATION_FREQUENCY:
	    default:
		max = (1 - sumCoalescentContributionRecombination * transitionProbabilityParameters.get(TransitionProbabilityParameters.ParameterChoice.RECOMBINATION_FREQUENCY)) / sumCoalescentContributionHybridization;
		break;
	    }

	    // update appropriately
	    if (!initializedFlag || (max < minmax)) {
		minmax = max;
		initializedFlag = true;
	    }
	}

	return (minmax);
    }

    /**
     * By construction, rows of a_ij matrix sum to one.
     */
    protected boolean verifyAij (double[][] a) {
	for (int i = 0; i < a.length; i++) {
	    for (int j = 0; j < a[i].length; j++) {
		if ((a[i][j] < 0.0) || (a[i][j] > 1.0)) {
		    System.err.println ("ERROR: entry in a_ij transition matrix is invalid. " + a[i][j]);
		    return (false);
		}
	    }
	}

	return (true);
    }

    /**
     * Call this after changes to parental tree branches, which affect
     * coalescent model calculations and the transition probabilities.
     *
     * Hmm... shall we do E-M learning in this class?
     */
    public void updateTransitionProbabilities () {
	double[] pi = calculatePi();
	double[][] a = calculateAij();
	
	for (int i = 0; i < pi.length; i++) {
	    myhmm.setPi(i, pi[i]);
	}

	for (int i = 0; i < a.length; i++) {
	    for (int j = 0; j < a[i].length; j++) {
		myhmm.setAij(i, j, a[i][j]);
	    }
	}
    }

    /**
     * Build HMM
     */
    protected void buildInitialHMM () {
	System.out.println("\n\nNow building initialHMM . . .");
	double[] pi = calculatePi();
	double[][] a = calculateAij();

	// kliu testing
	System.out.println ("Initial pi values:");
	for (double piv : pi) {
	    System.out.print (piv + " ");
	}
	System.out.println();
	
	System.out.println ("Initial a_ij values: ");
	for (int i = 0; i < a.length; i++) {
	    for (int j = 0; j < a[i].length; j++) {
		System.out.print (a[i][j] + " " );
	    }
	    System.out.println();
	}
	System.out.println();

	myhmm = buildMyHmm(trees_states, pi, a);
	System.out.println(myhmm);
    }

//(int) Math.pow(fParser.getAlphabet().size(),fParser.getNumSeq()
	
    /**
     * Read and parse basic info file
     * to build the Parser for reading sequences
     * @throws Exception 
     */
    protected void buildParser() throws Exception {
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
    protected void buildAlleleSpeciesMap() throws Exception {
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
     * Convert an EvoTree into a PhyloNet Network object
     */
    protected Network<Double> convertPHMMTreeToPhyloNetNetwork (EvoTree etree) {
	ExNewickReader<Double> enr = new ExNewickReader<Double>(new StringReader(etree.toNewickString(true, true)));
	// kliu - change this over to a String data member
	// don't really use it anyways
	try {
	    Network<Double> network = enr.readNetwork();
	    return (network);
	}
	catch(IOException ioe) {
	    System.err.println(ioe);
	    ioe.printStackTrace();
	    return (null);
	}
    }

    /**
     * Convert an EvoTree into a PhyloNet Tree object.
     */
    protected Tree convertPHMMTreeToPhyloNetTree (EvoTree etree) {
	NewickReader nr = new NewickReader(new StringReader(etree.toNewickString(true, true)));
	// kliu - change this over to a String data member
	// don't really use it anyways
	STITree<Double> newtr = new STITree<Double>(true);
	try {
	    nr.readTree(newtr);
	    return (newtr);
	}
	catch(Exception e) {
	    System.err.println(e);
	    e.printStackTrace();
	    return (null);
	}
    }
    	
    /**
     * Read and parse trees file
     * @throws Exception 
     */
    protected void buildTrees() throws Exception {
	if (parentalTreesFileName == null) {
	    throw new ParserFileException("Cannot read Trees file!");
	}

	trees_states = new ArrayList<HiddenState>();
	// also maintain equivalence classes among hidden states based on shared parental trees
	parentalTreeClasses = new HashMap<Network<Double>,Set<HiddenState>>();
	// also retain parental tree names
	// to facilitate parameter inputs/constraints on parental tree branches
	parentalTreeNameMap = new BijectiveHashtable<String,Network<Double>>();

	System.out.println("\nNow building trees . . .");
	BufferedReader ptreesbr = new BufferedReader(new FileReader(parentalTreesFileName));
	TreeParser ptp = new TreeParser(ptreesbr);
	ArrayList<EvoTree> eParentalTrees = ptp.nexusFileTreeNames(parentalTreesFileName);
	ptreesbr.close();
	Vector<Network<Double>> parentalTrees = new Vector<Network<Double>>();
	for (EvoTree etree : eParentalTrees) {
	    Network<Double> parentalTree = convertPHMMTreeToPhyloNetNetwork(etree);
	    // no duplicate node names allowed!
	    if (parentalTree.hasDuplicateNames()) {
		throw (new RuntimeException("ERROR: duplicate node names are present in parental tree " + etree.getName() + ". Check inputs and try again."));
	    }
	    parentalTrees.add(parentalTree);
	    // no duplicate parental tree names allowed!
	    if (parentalTreeNameMap.containsKey(etree.getName())) {
		throw (new RuntimeException("ERROR: duplicate parental tree name " + etree.getName() + ". Check inputs and try again."));
	    }
	    parentalTreeNameMap.put(etree.getName(), parentalTree);
	}
	
	// kliu - indexing is by (parentalTree, geneGenealogy) appearance order according to the following:
	for (Network<Double> parentalTree : parentalTrees) {
	    // kliu - cheap hack to clone all gene genealogies across
	    // each parental tree
	    //
	    // for now, don't share gene genealogies across parental trees (e.g. branch lengths)
	    BufferedReader ggbr = new BufferedReader(new FileReader(geneGenealogiesFileName));
	    TreeParser gtp = new TreeParser(ggbr);
	    ArrayList<EvoTree> eGeneGenealogies = gtp.nexusFileTreeNames(geneGenealogiesFileName);
	    ggbr.close();

	    HashSet<HiddenState> parentalTreeEquivalenceClass = new HashSet<HiddenState>();
	    // disallow duplicate gene genealogy names
	    HashSet<String> geneGenealogyNames = new HashSet<String>();

	    for (EvoTree egg : eGeneGenealogies) {
		Tree geneGenealogy = convertPHMMTreeToPhyloNetTree(egg);
		if (hasDuplicateNames(geneGenealogy)) {
		    throw (new RuntimeException("ERROR: duplicate node names are present in gene genealogy " + egg.getName() + ".Check inputs and try again."));
		}
		String hiddenStateName = parentalTreeNameMap.rget(parentalTree) + HiddenState.HIDDEN_STATE_NAME_DELIMITER + egg.getName();
		// kliu - meh - parse allele-to-species mapping later and add in references here
		HiddenState hiddenState = new HiddenState(hiddenStateName, parentalTree, geneGenealogy, null, parentalTreeEquivalenceClass);
		trees_states.add(hiddenState);
		parentalTreeEquivalenceClass.add(hiddenState);
		// really strict
		if (parentalTreeNameMap.containsKey(egg.getName())) {
		    throw (new RuntimeException("ERROR: gene genealogy name " + egg.getName() + " appears as a parental tree name. Check inputs and try again."));
		}

		if (geneGenealogyNames.contains(egg.getName())) {
		    throw (new RuntimeException("ERROR: duplicate gene genealogy name " + egg.getName() + ". Check inputs and try again."));
		}
		geneGenealogyNames.add(egg.getName());
	    }

	    // maintain the map
	    parentalTreeClasses.put(parentalTree, parentalTreeEquivalenceClass);
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

	// strict!
	if (trees_states.size() < 2) {
	    throw (new ParserFileException("ERROR: must be at least two hidden states in PhyloNet-HMM."));
	}
    }
    
    /**
     * Read transition probability parameters (other than those related to basic coalescent model calculations, i.e., parental tree
     * branch lengths).
     *
     * Current model uses two parameters: a hybridization frequency $v$, and a recombination frequency $u$.
     * See writeup for details.
     */
    protected void readTransitionProbabilityParameters (BufferedReader br) throws Exception {
	System.out.println("Input transition probability parameters in format <recombination frequency parameter $u$> <hybridization frequency parameter $v$>");
	String line = br.readLine();
	StringTokenizer st = new StringTokenizer(line);
	if (st.countTokens() != 2) {
	    throw (new IOException("ERROR: invalid transition probability parameters."));
	}
	double recombinationFrequency = Double.parseDouble(st.nextToken());
	double hybridizationFrequency = Double.parseDouble(st.nextToken());
	transitionProbabilityParameters = new TransitionProbabilityParameters(recombinationFrequency, hybridizationFrequency);
    }
}










    // /**
    //  * Get Initial Pi Probabilities array
    //  * @throws IOException - rare
    //  */
    // protected static void getPiInfo(BufferedReader in) throws IOException {
    // 	boolean getPi = true;
		
    // 	while (getPi) {
			
    // 	    // Getting pi
    // 	    System.out.println("\nThe Pi or Initial Probabilities. \n(Note: Make sure pi probabilities correspond to the same order of trees in the trees input file.) \n Sample input for 4 trees/states: .2 .3 .1 .4");
    // 	    System.out.println("Input Initial Pi probabilities:");
    // 	    String[] piString;
    // 	    piString = in.readLine().split(" ");

				
    // 	    // Error Check: check to see if number of states and number of pi probabilities are the same
    // 	    if (numStates == piString.length) {
    // 		getPi = false;
					
    // 		// make pi probabilities array
    // 		pi = new double[numStates];
    // 		double sum = 0;
    // 		for (int i=0; i < piString.length; i++) {
    // 		    try {
    // 			pi[i] = Double.parseDouble(piString[i]);
    // 			sum += pi[i];
    // 		    } 
    // 		    catch (NumberFormatException e) { 
    // 			getPi = true;
    // 			System.out.println("Number format error: Cannot convert inputed pi number : " + piString[i] + " to a double."); 
    // 		    }
    // 		}
    // 		if (Math.abs(1.0 - sum) > tolerated_error) {
    // 		    System.out.println("Error: the sum of the pi probabilities did not sum up to 1.0");
    // 		    getPi = true;
    // 		}
					
    // 	    } else System.out.println("Error: number of pi probability inputs did not match number of states! ");

    // 	}
    // }




    // /**
    //  * Get State Transition Matrix Aij
    //  * @throws IOException - rare
    //  */
    // protected static void getAij(BufferedReader in) throws IOException {
    // 	boolean getA = true;
    // 	boolean getChoice = true;
		
    // 	while (getA) {			
    // 	    // Getting transition matrix A
    // 	    System.out.println("\nThe Transition of States/Trees Matrix Aij. \n(Note: Make sure the order in the matrix correspond to the same order of trees in the given tree input file.");
    // 	    while (getChoice) {
    // 		System.out.println("Choose an option: \n(a) Read transition matrix by file.\n(b) Input transition matrix manually.");
    // 		String choice = in.readLine().toLowerCase();
    // 		if (choice.equals("a")) {
    // 		    try {
    // 			System.out.println("Please see README for .matrix file format.");
    // 			System.out.println("Input .matrix file path name:");
    // 			String filename = in.readLine();
    // 			BufferedReader filebr = new BufferedReader(new FileReader(filename));
    // 			aij = readAij(numStates, filebr);
    // 			filebr.close();
    // 			getChoice = false;
    // 		    } catch (Exception e) {
    // 			System.out.println(e);
    // 			System.out.println("Please try again. \n");
    // 		    }
    // 		}
    // 		else if (choice.equals("b")) {
    // 		    try {
    // 			System.out.println("Example transition for 3 trees (inputted 3 times): \n   .3 .2 .5 \n   .2 .4 .4 \n   .5 .4 .1 \n");
    // 			System.out.println("Input transition Aij matrix:");
    // 			aij = readAij(numStates, in);
    // 			getChoice = false;
    // 		    } catch (Exception e) {
    // 			System.out.println(e);
    // 			System.out.println("Please try again.");
    // 		    }
    // 		}
    // 	    }
    // 	    getA = false;
    // 	}
    // }
	
	
	
	
	
	
    // /**
    //  * Read in Transition matrix
    //  * @param numStates - the number of states/trees
    //  * @param br - a buffered reader
    //  * @returns a double[][] transition of states/trees matrix probabilities
    //  * @throws Exception - Either an IO exception if an IO error occurred or a ParserFileException that signals an incorrect matrix input.
    //  */
    // public static double[][] readAij(int numStates, BufferedReader br) throws Exception {
    // 	String[] row;
    // 	double[] newRow;
    // 	double[][] aij = new double[numStates][numStates];
		
    // 	for (int i = 0; i < numStates; i++) {
    // 	    newRow = new double[numStates];
    // 	    row = br.readLine().split(" ");
    // 	    double sum = 0;
    // 	    if (row.length != numStates) throw new ParserFileException("Error: Number of entries in matrix is incorrect!");
    // 	    for (int j = 0; j < row.length; j++) {
    // 		newRow[j] = Double.parseDouble(row[j]);
    // 		sum += newRow[j];
    // 	    }
			
    // 	    // Check to make sure the sum of the entire row sums up to 1.0
    // 	    if (Math.abs(1.0 - sum) > tolerated_error) {
    // 		throw new ParserFileException("Error: the sum of row " + i + " did not sum up to 1.0.");
    // 	    }
			
    // 	    aij[i] = newRow;
    // 	}
		
    // 	// Check to make sure the sum of each column sums up to 1.0
    // 	for (int j = 0; j < numStates; j++) {
    // 	    double sum = 0;
    // 	    for (int i = 0; i < numStates; i++) {
    // 		sum+= aij[i][j];
    // 	    }
    // 	    if (Math.abs(1.0 - sum) > tolerated_error) {
    // 		throw new ParserFileException("Error: the sum of col " + j + " did not sum up to 1.0.");
    // 	    }
    // 	}
    // 	return aij;
    // }
