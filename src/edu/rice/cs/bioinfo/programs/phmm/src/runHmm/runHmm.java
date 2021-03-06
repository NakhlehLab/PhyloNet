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

package edu.rice.cs.bioinfo.programs.phmm.src.runHmm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.StringWriter;
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
import java.util.Collections;

import edu.rice.cs.bioinfo.programs.phmm.src.util.Constants;
import edu.rice.cs.bioinfo.programs.phmm.src.util.Matrix;
import edu.rice.cs.bioinfo.programs.phmm.src.util.TreeUtils;
import edu.rice.cs.bioinfo.programs.phmm.src.util.MapOfMap;
import edu.rice.cs.bioinfo.programs.phmm.src.substitutionModel.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phmm.src.substitutionModel.GTRSubstitutionModel;
import edu.rice.cs.bioinfo.programs.phmm.src.substitutionModel.NucleotideAlphabet;
import edu.rice.cs.bioinfo.programs.phmm.src.optimize.MultivariateOptimizer;
import edu.rice.cs.bioinfo.programs.phmm.src.optimize.CalculationCache;
//import gridSearch.GridSearchAlgorithm;
import edu.rice.cs.bioinfo.programs.phmm.src.phylogeny.EvoTree;
import edu.rice.cs.bioinfo.programs.phmm.src.phylogeny.TreeParser;
import edu.rice.cs.bioinfo.programs.phmm.src.reader.Parser;
import edu.rice.cs.bioinfo.programs.phmm.src.reader.ParserFileException;
import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.Hmm;
import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm.OpdfMap;
import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm.SwitchingFrequency;
import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm.SwitchingFrequencyRatioTerm;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF.CoalescePattern;
import edu.rice.cs.bioinfo.library.programming.BijectiveHashtable;
import edu.rice.cs.bioinfo.library.programming.BidirectionalMultimap;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;

public class runHmm {
    protected static final String SPECIES_TO_ALLELES_MAPPING_ENTRY_DELIMITER = ",";
    protected static final String SPECIES_TO_ALLELES_MAPPING_KEY_VALUE_DELIMITER = ":";
    protected static final String SWITCHING_FREQUENCY_RATIO_TERM_TRANSITION_NAMES_DELIMITER_1 = "|";
    protected static final String SWITCHING_FREQUENCY_RATIO_TERM_TRANSITION_NAMES_DELIMITER_2 = ",";

    protected static final double tolerated_error = 1e-5;		/* Sum of probabilities margin of error allowed */
	
    // stored information 
    protected String basicFileName = null;						/* Basic Info file name */
    protected String parentalTreesFileName = null;				/* Parental trees file name */
    protected String geneGenealogiesFileName = null;			/* Gene genealogies file name */
    
    // kliu - also need an empty temporary working directory
    protected String workingDirectory;

    
    // information created/built
    // kliu - index into tuples
    protected ArrayList<HiddenState> hiddenStates;			/* List of all states/trees */
    protected Hmm<ObservationMap> myhmm;			/* The entire HMM */

    // argh - in lieu of worrying about EvoTree.equals() method
    // Maintain equivalence classes among hidden states based on shared parental tree.
    protected Map<Network<CoalescePattern[]>,Set<HiddenState>> parentalTreeClasses;
    // ditto for gene genealogy
    protected Map<Tree,Set<HiddenState>> rootedGeneGenealogyClasses;
    protected Map<Tree,Set<HiddenState>> unrootedGeneGenealogyClasses;
    // Since it's bijective, make it reversible
    protected BijectiveHashtable<String,Network<CoalescePattern[]>> parentalTreeNameMap;
    // bleh, need to keep track of allele-to-species mapping for
    // each parental tree
    //
    // Hmm... Map.equals() function may be too clever here.
    // Just use a Map.
    protected Map<String,Map<String,List<String>>> parentalTreeSpeciesToAllelesMapMap;
    // ditto for gene genealogies
    protected BijectiveHashtable<String,Tree> geneGenealogyNameMap;
    // only topologies matter for rootedGeneGenealogies,
    // and unrooted representative of a topological equivalence class
    // has its branch lengths optimized.
    protected BidirectionalMultimap<Tree,Tree> rootedToUnrootedGeneGenealogyMap;
    // name of outgroup taxon
    protected String outgroupTaxonName;

    // similarly, need a map from (parental tree name, gene genealogy name) -> hidden state object
    // since need to specify in hidden state switching frequency ratio file
    protected MapOfMap<String,String,HiddenState> parentalTreeGeneGenealogyNamePairToHiddenStateMap;

    //protected boolean collapseGeneGenealogiesByTopologicalEquivalenceClassFlag;

    // External state to map
    // (Gene genealogy name 1, gene genealogy name 2)->SwitchingFrequencyRatioTerm object for gene genealogy pairs.
    protected BidirectionalMultimap<Tuple<HiddenState,HiddenState>,SwitchingFrequencyRatioTerm> hiddenStatePairToSwitchingFrequencyRatioTermMap;
    protected BijectiveHashtable<String,SwitchingFrequencyRatioTerm> nameToSwitchingFrequencyRatioTermMap;
    protected Hashtable<String,Tuple3<Double,Double,Boolean>> switchingFrequencyRatioTermNameToOptimizeFlagMap;

    protected SwitchingFrequency gamma;
    
    // need to assume GTR substitution model due to parameterization differences
    protected GTRSubstitutionModel gtrSubstitutionModel;

    // single cache shared amongst all objects
    protected CalculationCache calculationCache;

    // kliu - neither are necessary anymore
    //    public static double[] pi;								/* The initial pi probabilities for each state */    
    //    public static double[][] aij;							/* The state transition probability matrix */
    protected Parser fParser;							/* The parser for all basic info and read sequences --> also calculates likelihood */
    protected int numStates = -1;						/* The number of states for the HMM */
	
    protected static void printUsage () {
	System.err.println ("Prompt-based usage: java -jar dist/lib/phmm.jar");
	System.err.println ("File-based usage: java -jar dist/lib/phmm.jar <text file with input commands>");
    }

    public runHmm () {
	calculationCache = new CalculationCache();
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

    public String getWorkingDirectory () {
	return (workingDirectory);
    }

    /**
     * WARNING: returns null to indicate no outgroup taxon available.
     */
    public String getOutgroupTaxonName () {
	return (outgroupTaxonName);
    }

    /**
     * Since some construction performed by MultivariateOptimizer.
     */
    // public boolean getCollapseGeneGenealogiesByTopologicalEquivalenceClassFlag () {
    // 	return (collapseGeneGenealogiesByTopologicalEquivalenceClassFlag);
    // }

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
	    
	    // System.out.println ("Collapse gene genealogies according to topological equivalence classes? Enter true for yes or anything else for no:" );
	    // collapseGeneGenealogiesByTopologicalEquivalenceClassFlag = Boolean.parseBoolean(in.readLine());
			
	    System.out.println("\nInput the parental trees file path name:\n (note: see README for file format) \n");
	    parentalTreesFileName = in.readLine();
			
	    System.out.println("\nInput the gene genealogies file path name:\n (note: see README for file format) \n");
	    geneGenealogiesFileName = in.readLine();

	    System.out.println("\nInput outgroup taxon name, or empty string for no outgroup taxon: \n");
	    outgroupTaxonName = in.readLine().trim();
	    if (outgroupTaxonName.equals("")) {
		outgroupTaxonName = null;
	    }
	    
	    System.out.println("Empty working directory: ");
	    workingDirectory = in.readLine();

	    readSubstitutionModelParameters(in);

	    // Get the Pi probabilities array
	    //getPiInfo(in);
	    
	    // Get transition probability parameters.
	    // See writeup for details.
	    //getAij(in);
			
	    // Now build the trees and get the tree mapping to integers
	    buildTrees();
	    
	    //Reading in Basic Info file and store information
	    buildParser();
		
	    // Bleh - need to read in and construct TransitionProbabilityParameter objects here.
	    // Specified according to parental-tree-switching-parameter and gene-genealogy-switching-parameter
	    // file discussed in README.
	    // Then, in MultivariateOptimizer, add SwitchingFrequencyParameter for each perform actual optimization
	    // for each TransitionProbabilityParameter.
	    // add map to parameter
	    readSwitchingFrequency(in);

	    // need to push this to after the trees are read in
	    // since term maximums depend on the number of 
	    // parental/gene trees
	    readSwitchingFrequencyRatioTermFile(in);

	    // verify switching frequency ratio terms
	    if (!verifySwitchingFrequencyRatioTerms()) {
		System.err.println ("ERROR: unable to verify switching frequency ratio terms. Please correct and try again.");
		System.exit(1);
	    }

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
		// ArrayList<ObservationMap> obsSeq;
		// int seqChoice = sequenceChoice(in);
		// if (seqChoice == 1) {
		//     obsSeq = getObs(in);
		// } else {
		//     obsSeq = fParser.getObs();
		// }
		
		ArrayList<ArrayList<ObservationMap>> obsSeqArray = getDataPartitions(in);

		Vector<Double> partitionViterbiLikelihoods = new Vector<Double>();
		Vector<Double> partitionModelLikelihoods = new Vector<Double>();

		// log likelihoods
		double viterbiLLH = 0.0;
		double llh = 0.0;

		for (int partition = 0; partition < obsSeqArray.size(); partition++) {
		    Tuple<int[],Double> viterbiResult = myhmm.viterbiStateSequence(obsSeqArray.get(partition));	
		    double partitionModelLikelihood = myhmm.lnProbability(obsSeqArray.get(partition));

		    BufferedWriter bw = new BufferedWriter(new FileWriter(outputfile + MultivariateOptimizer.FILENAME_SUFFIX_DELIMITER + Integer.toString(partition)));
		    for (int index : viterbiResult.Item1) {
			bw.write(hiddenStates.get(index).getName()); bw.newLine();
		    }
		    bw.flush();
		    bw.close();

		    viterbiLLH += viterbiResult.Item2.doubleValue();
		    llh += partitionModelLikelihood;
		}

		System.out.println ("Input HMM Viterbi log likelihood: |" + viterbiLLH + "|");
		System.out.println ("Input HMM log likelihood: |" + llh + "|");
	
		break;
	    // case 1:
	    // 	// RUN BAUM WELCH
	    //     // Modifies the current Hmm;

	    //     System.out.println("\n");

	    // 	// Testing purposes
	    // 	//			ArrayList<ObservationMap> oSeq = listOfOseq.get(0);
	    // 	//			System.out.println("\n ---------------------\nThe HMM Emission Probabilities");
	    // 	//			for (int i = 0; i < hiddenStates.size(); i++) {
	    // 	//				System.out.println("State : " + i);
	    // 	//				System.out.println("Opdf : {");
	    // 	//				for (int j = 0; j < oSeq.size(); j++) {
	    // 	//					System.out.print(oSeq.get(j) + ":" + myhmm.getOpdf(i).probability(oSeq.get(j)) + " , ");
	    // 	//				}
	    // 	//				System.out.print("}\n\n");
	    // 	//			}
			
	    // 	// END TEST PURPOSES

	    //     // Allow users to read in New Sequence Observation or re-use a
	    //     // previously read in Sequence
	    //     int seqChoice2 = sequenceChoice(in);
	    //     ArrayList<ObservationMap> obsSequence;
	    //     if (seqChoice2 == 1) {
	    //         obsSequence = getObs(in);
	    //     } else {
	    //         obsSequence = fParser.getObs();
	    //     }

	    // 	BaumWelchScaledLearner bwl = new BaumWelchScaledLearner();

	    // 	Hmm<ObservationMap> learnedhmm = bwl.learn(myhmm, obsSequence, 9);
	    // 	myhmm = learnedhmm;
	    // 	learnedhmm = null; // for garbage collection
	    // 	System.out.println("\n----------------------------\nThe LEARNED HMM IS: ");
	    // 	System.out.println(myhmm + "\n\n");
			
	    // 	break;
		
		
	    // case 2:
	    //     // RUN GRID SEARCH
	    //     // Modifies the current hmm;

	    //     System.out.println("\n");
	    //     // Allow users to read in New Sequence Observation or re-use a
	    // 	// previously read in Sequence
	    // 	int seqChoice3 = sequenceChoice(in);
	    // 	ArrayList<ObservationMap> seqObs;
	    // 	if (seqChoice3 == 1) {
	    // 	    seqObs = getObs(in);
	    // 	} else {
	    // 	    seqObs = fParser.getObs();
	    // 	}

	    // 	//Get the max and mins and g samples of each parameter
	    // 	try {
	    // 	    System.out.println("Input the number of samples for tree branch lengths: ");
	    // 	    int gBranchIn = Integer.parseInt(in.readLine());

	    // 	    System.out.println("Input the number of samples for Recombination Frequency: ");
	    // 	    int gRecombinationIn = Integer.parseInt(in.readLine());

	    // 	    System.out.println("Input the number of samples for Hybridization Frequency: ");
	    // 	    int gHybridizationIn = Integer.parseInt(in.readLine());

	    // 	    System.out.println("Input the number of samples for Felsenstein Base Substitution Rate: ");
	    // 	    int gBaseSubIn = Integer.parseInt(in.readLine());

	    // 	    System.out.println("Input the minimum double value for all tree branch lengths: ");
	    // 	    double branchMinIn = Double.parseDouble(in.readLine());

	    // 	    System.out.println("Input the maximum double value for all tree branch lengths: ");
	    // 	    double branchMaxIn = Double.parseDouble(in.readLine());

	    // 	    // System.out.println("Input the minimum double value for Recombination frequency: ");
	    // 	    // double recombinationMinIn = Double.parseDouble(in.readLine());

	    // 	    // System.out.println("Input the maximum double value for Recombination frequency: ");
	    // 	    // double recombinationMaxIn = Double.parseDouble(in.readLine());

	    // 	    System.out.println("Input the minimum double value for Hybridization frequency: ");
	    // 	    double hybridizationMinIn = Double.parseDouble(in.readLine());

	    // 	    System.out.println("Input the maximum double value for Hybridization frequency: ");
	    // 	    double hybridizationMaxIn = Double.parseDouble(in.readLine());

	    // 	    System.out.println("Input the minimum double value for Felsenstein base substitution rate: ");
	    // 	    double baseSubMinIn = Double.parseDouble(in.readLine());

	    // 	    System.out.println("Input the maximum double value for Felsenstein base substitution rate: ");
	    // 	    double baseSubMaxIn = Double.parseDouble(in.readLine());

	    // 	    // gRecombinationIn,
	    // 	    // recombinationMinIn, recombinationMaxIn, 
	    // 	    GridSearchAlgorithm gsa =
            //             new GridSearchAlgorithm(gBranchIn, 
	    // 					gHybridizationIn, gBaseSubIn, branchMinIn, branchMaxIn,
	    // 					hybridizationMinIn,
	    // 					hybridizationMaxIn, baseSubMinIn, baseSubMaxIn, this);

	    // 	    gsa.runGridSearch(seqObs, myhmm, transitionProbabilityParameters,
	    // 			      hiddenStates, parentalTreeClasses);

	    // 	}
	    // 	catch (Exception e) {
	    // 	    System.err.println("Input error : " + e);
	    // 	    break;
	    // 	}

	    //     break;
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
	    // int sc = sequenceChoice(in);
	    // ArrayList<ObservationMap> obsSequence;
	    // if (sc == 1) {
	    // 	obsSequence = getObs(in);
	    // } getData
	    // else {
	    // 	obsSequence = fParser.getObs();
	    // }

	    // meh, just read in new sequence for this operation
	    ArrayList<ArrayList<ObservationMap>> obsSequenceArray = getDataPartitions(in);	
	    
	    System.out.println("Input parental-branch-length-parameter-to-edge map filename: ");
	    String inputLengthParameterToEdgeMapFilename = in.readLine();
	    System.out.println("Input parental-branch-length-parameter strict inequalities filename: ");
	    String inputParentalBranchLengthParameterInequalitiesFilename = in.readLine();
	    System.out.println("Input length-parameter-set-constraints filename: ");
	    String inputLengthParameterSetConstraintsFilename = in.readLine();
	    System.out.println("Input checkpoint file to restore from, or empty line for no restore: ");
	    String inputRestoreCheckpointFilename = in.readLine();
	    System.out.println("Output posterior decoding probabilities file: ");
	    String posteriorDecodingProbabilitiesFilename = in.readLine();
	    System.out.println("Output Viterbi-optimal hidden state sequence file: ");
	    String viterbiHiddenStateSequenceFilename = in.readLine();
	    System.out.println("Output model likelihoods file: ");
	    String modelLikelihoodsFilename = in.readLine();
	    System.out.println("Output file with optimized model parameter values: ");
	    String optimizedModelParameterValuesFilename = in.readLine();
	    System.out.println("Initial search settings vector <setting 1> <setting 2> ... <setting s>, where <setting s> is one of CURRENT, DEFAULT, RANDOM");
	    MultivariateOptimizer.InitialSearchSettings[] initialSearchSettings = parseInitialSearchSettings(in.readLine());
	    System.out.println("Enable optimization flag vector <enable parental tree optimization flag> <enable gene genealogy optimization flag> <enable switching frequency optimization flag> <enable substitution model optimization flag>");
	    StringTokenizer flagTokens = new StringTokenizer(in.readLine());
	    if (flagTokens.countTokens() != 4) { throw (new RuntimeException("ERROR: invalid optimization flag vector.")); }
	    boolean enableParentalTreeOptimizationFlag = Boolean.parseBoolean(flagTokens.nextToken());
	    boolean enableGeneGenealogyOptimizationFlag = Boolean.parseBoolean(flagTokens.nextToken());
	    boolean enableSwitchingFrequencyOptimizationFlag = Boolean.parseBoolean(flagTokens.nextToken());
	    boolean enableSubstitutionModelOptimizationFlag = Boolean.parseBoolean(flagTokens.nextToken());

	    MultivariateOptimizer multivariateOptimizer = new MultivariateOptimizer(myhmm,
										    this,
										    hiddenStates,
										    gtrSubstitutionModel,
										    parentalTreeNameMap,
										    geneGenealogyNameMap,
										    rootedToUnrootedGeneGenealogyMap,
										    gamma,
										    nameToSwitchingFrequencyRatioTermMap,
										    switchingFrequencyRatioTermNameToOptimizeFlagMap,
										    obsSequenceArray,
										    inputLengthParameterToEdgeMapFilename,
										    inputParentalBranchLengthParameterInequalitiesFilename,
										    inputLengthParameterSetConstraintsFilename,
										    calculationCache,
										    enableParentalTreeOptimizationFlag,
										    enableGeneGenealogyOptimizationFlag,
										    enableSwitchingFrequencyOptimizationFlag,
										    enableSubstitutionModelOptimizationFlag,
										    inputRestoreCheckpointFilename
										    );

	    System.out.println ("Optimizing PhyloNet-HMM parameters... ");
	    multivariateOptimizer.optimize(initialSearchSettings);
	    System.out.println ("Optimizing PhyloNet-HMM parameters DONE. ");

	    System.out.println ("Saving posterior decoding probabilities... "); 
	    for (int partition = 0; partition < obsSequenceArray.size(); partition++) {
		BufferedWriter pdpbw = new BufferedWriter(new FileWriter(posteriorDecodingProbabilitiesFilename + MultivariateOptimizer.FILENAME_SUFFIX_DELIMITER + Integer.toString(partition)));
		// indexed by [observation index][hidden state index]
		double[][] posteriorDecodingProbabilities = multivariateOptimizer.computeHMMPosteriorDecodingProbabilities(obsSequenceArray.get(partition));
		for (int i = 0; i < posteriorDecodingProbabilities.length; i++) {
		    for (int j = 0; j < posteriorDecodingProbabilities[i].length; j++) {
			pdpbw.write (i + " " + j + " " + posteriorDecodingProbabilities[i][j]); pdpbw.newLine();
		    }
		}
		pdpbw.flush();
		pdpbw.close();
	    }
	    System.out.println ("Saving posterior decoding probabilities DONE."); 

	    double viterbiLogLikelihood = 0.0;
	    double modelLogLikelihood = 0.0;
	    System.out.println ("Saving Viterbi-optimal hidden state sequence and computing Viterbi-optimal and model log likelihoods... "); 
	    for (int partition = 0; partition < obsSequenceArray.size(); partition++) {
		Tuple<int[],Double> viterbiResult = myhmm.viterbiStateSequence(obsSequenceArray.get(partition));
		BufferedWriter bw = new BufferedWriter(new FileWriter(viterbiHiddenStateSequenceFilename + MultivariateOptimizer.FILENAME_SUFFIX_DELIMITER + Integer.toString(partition)));
		viterbiLogLikelihood += viterbiResult.Item2;
		for (int index : viterbiResult.Item1) {
		    bw.write(hiddenStates.get(index).getName()); bw.newLine();
		}
		bw.flush();
		bw.close();

		modelLogLikelihood += myhmm.lnProbability(obsSequenceArray.get(partition));
	    }
	    System.out.println ("Saving Viterbi-optimal hidden state sequence and computing Viterbi-optimal and model log likelihoods DONE."); 

	    System.out.println ("Optimized model's final log likelihood: |" + modelLogLikelihood + "|");
	    System.out.println ("Optimized model's Viterbi-optimal log likelihood: |" + viterbiLogLikelihood + "|");

	    System.out.println ("Saving log likelihoods... ");
	    BufferedWriter lbw = new BufferedWriter(new FileWriter(modelLikelihoodsFilename));
	    lbw.write(Double.toString(modelLogLikelihood)); lbw.newLine();
	    lbw.write(Double.toString(viterbiLogLikelihood)); lbw.newLine();
	    lbw.flush();
	    lbw.close();
	    System.out.println ("Saving log likelihoods DONE. ");
	    
	    // also save values for all model parameters
	    outputOptimizedModelParameterValues(optimizedModelParameterValuesFilename);
	}
	catch (IOException ioe) {
	    System.err.println (ioe);
	    ioe.printStackTrace();
	}
    }

    protected MultivariateOptimizer.InitialSearchSettings[] parseInitialSearchSettings (String line) {
	StringTokenizer st = new StringTokenizer(line);
	Vector<MultivariateOptimizer.InitialSearchSettings> vec = new Vector<MultivariateOptimizer.InitialSearchSettings>();
	while (st.hasMoreTokens()) {
	    vec.add(MultivariateOptimizer.InitialSearchSettings.valueOf(st.nextToken()));
	}
	MultivariateOptimizer.InitialSearchSettings[] result = new MultivariateOptimizer.InitialSearchSettings[vec.size()];
	result = vec.toArray(result);
	return (result);
    }

    protected void outputOptimizedModelParameterValues (String filename) {
	try {
	    // bleh
	    RnNewickPrinter<CoalescePattern[]> rnNewickPrinter = new RnNewickPrinter<CoalescePattern[]>();
	    BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
	    bw.write("Hidden state order: "); bw.newLine();
	    for (int i = 0 ; i < hiddenStates.size(); i++) {
		bw.write(i + " " + hiddenStates.get(i).getName()); bw.newLine();
	    }
	    bw.newLine();
	    List<String> parentalTreeNames = new ArrayList<String>(parentalTreeNameMap.keys());
	    Collections.sort(parentalTreeNames);
	    for (String parentalTreeName : parentalTreeNames) {
		bw.write("Parental tree " + parentalTreeName + ": "); bw.newLine();
		// bleh again
		StringWriter sw = new StringWriter();
		rnNewickPrinter.print(parentalTreeNameMap.get(parentalTreeName), sw);
		bw.write(sw.toString()); bw.newLine();
		bw.newLine();
	    }
	    bw.newLine();
	    for (HiddenState hiddenState : hiddenStates) {
		bw.write("Rooted gene genealogy associated with hidden state " + hiddenState.getName() + ":"); bw.newLine();
		bw.write(hiddenState.getRootedGeneGenealogy().toNewick()); bw.newLine();
		bw.newLine();
	    }
	    bw.newLine();
	    for (HiddenState hiddenState : hiddenStates) {
		bw.write("Unrooted gene genealogy associated with hidden state " + hiddenState.getName() + ":"); bw.newLine();
		bw.write(hiddenState.getUnrootedGeneGenealogy().toNewick()); bw.newLine();
		bw.newLine();
	    }
	    bw.newLine();
	    for (HiddenState hiddenState : hiddenStates) {
		bw.write("Processed rooted gene genealogy associated with hidden state " + hiddenState.getName() + ":"); bw.newLine();
		bw.write(hiddenState.getProcessedRootedGeneGenealogy().toNewick()); bw.newLine();
		bw.newLine();
	    }
	    bw.newLine();
	    for (String name : nameToSwitchingFrequencyRatioTermMap.keys()) {
		bw.write("Hidden-state-switching frequency ratio term parameter with name " + name + ": " + nameToSwitchingFrequencyRatioTermMap.get(name).getValue()); bw.newLine();
	    }
	    bw.newLine();
	    bw.write("Hidden-state-switching frequency with name " + gamma.getName() + ": " + gamma.getValue()); bw.newLine();
	    bw.newLine();
	    bw.write("Normalized within-row hidden-state-switching frequencies f'_{ijk}: "); bw.newLine();
	    // for readability, also output corresponding switching frequencies
	    // just need to multiply through by (1-\gamma) or \gamma / (num-parental-trees - 1), as appropriate
	    MapOfMap<HiddenState,HiddenState,Double> hiddenStateSwitchingFrequencyMap = calculateSwitchingFrequencies();
	    //bw.write (hiddenStateSwitchingFrequencyMap.toString()); bw.newLine();
	    // for readability, proceed according to above hidden state order
	    for (HiddenState h1 : hiddenStates) {
		for (HiddenState h2 : hiddenStates) {
		    // later, if adding across-row switching frequency
		    // parameterization, remove this
		    double correctionFactor = gamma.getValue() / ((gamma.getNumAlternatives() - 1) * 1.0);
		    if (checkSameParentalTreeClass(h1, h2)) {
			correctionFactor = (1.0 - gamma.getValue());
		    }

		    // unnecessary by verifySwitchingFrequencyRatioTerms(),
		    // but paranoid
		    if (!hiddenStateSwitchingFrequencyMap.contains(h1, h2)) {
			System.err.println ("ERROR: in runHmm.outputOptimizedModelParameterValues(), missing transition for pair of hidden states named " + h1.getName() + " " + h2.getName() + ". Skipping.");
			continue;
		    }
		    bw.write(h1.getName() + " -> " + h2.getName() + " : " + (hiddenStateSwitchingFrequencyMap.get(h1, h2).doubleValue() / correctionFactor)); bw.newLine();
		}
	    }
	    bw.newLine();
	    
	    bw.write("GTR base frequencies: "); bw.newLine();
	    bw.write(Matrix.toString(gtrSubstitutionModel.getStationaryProbabilities())); bw.newLine();
	    bw.newLine();
	    bw.write("GTR substitution rate parameter values: "); bw.newLine();
	    bw.write(Matrix.toString(gtrSubstitutionModel.getOriginalRateParameters())); bw.newLine();
	    bw.newLine();
	    bw.write("GTR rate matrix: "); bw.newLine();
	    bw.write(Matrix.toString(gtrSubstitutionModel.getFullRateMatrix())); bw.newLine();
	    bw.newLine();
	    // bleh
	    double[][] a = new double[hiddenStates.size()][hiddenStates.size()];
	    for (int i = 0; i < hiddenStates.size(); i++) {
		for (int j = 0; j < hiddenStates.size(); j++) {
		    a[i][j] = myhmm.getAij(i, j);
		}
	    }
	    bw.write("PhyloNet-HMM transition probability matrix: "); bw.newLine();
	    bw.write(Matrix.toString(a)); bw.newLine();
	    bw.newLine();
	    bw.flush();
	    bw.close();
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
     * Set geneGenealogyTaxonFlag to true to get gene-genealogy taxa, or false to get parental-tree taxa.
     * Warning - constructs a new array each time.
     *
     * Gene-genealogy-taxon -> parental-tree-taxon mapping.
     */
    protected String[] getSortedTaxaFromSpeciesToAllelesMapping (Map<String,List<String>> speciesToAllelesMap, boolean geneGenealogyTaxonFlag) {
	// no official list of species names - meh
	Set<String> set;
	if (geneGenealogyTaxonFlag) {
	    set = new HashSet<String>();
	    for (String species : speciesToAllelesMap.keySet()) {
		for (String allele : speciesToAllelesMap.get(species)) {
		    set.add(allele);
		}
	    }
	}
	else {
	    set = new HashSet<String>(speciesToAllelesMap.keySet());
	    // special case for outgroup taxon
	    // sort of a hack
	    // meh
	    if ((getOutgroupTaxonName() != null) && (set.contains(getOutgroupTaxonName()))) {
		set.remove(getOutgroupTaxonName());
	    }
	}

	String[] sortedTaxa = new String[set.size()];
	sortedTaxa = set.toArray(sortedTaxa);
	Arrays.sort(sortedTaxa);
	return (sortedTaxa);
    }

    protected boolean verifyGeneGenealogyTaxa () {
	// unique values in species-alleles-mapping and in basic-info-file must be consistent with 
	// gene genealogies

	// basic info taxa
	String[] sortedBasicInfoGeneGenealogyTaxa = getSortedGeneGenealogyTaxa();

	// check (rooted) gene genealogy taxa against basic info taxa
	for (HiddenState hiddenState : hiddenStates) {
	    String[] sortedCheckTaxa = hiddenState.getRootedGeneGenealogy().getLeaves();
	    Arrays.sort(sortedCheckTaxa);
	    if (!Arrays.equals(sortedBasicInfoGeneGenealogyTaxa, sortedCheckTaxa)) {
		System.err.println ("ERROR: list of taxa in basic-info-file is not consistent with gene genealogies.");
		System.err.println (Arrays.toString(sortedBasicInfoGeneGenealogyTaxa) + "\n");
		System.err.println (Arrays.toString(sortedCheckTaxa) + "\n");
		return (false);
	    }
	}

	// check species-alleles-mapping against basic info taxa
	for (Map<String,List<String>> speciesToAllelesMap : parentalTreeSpeciesToAllelesMapMap.values()) {
	    for (HiddenState hiddenState : hiddenStates) {
		String[] sortedReferenceTaxa = getSortedTaxaFromSpeciesToAllelesMapping(hiddenState.getSpeciesToAllelesMapping(), true);
		String[] sortedCheckTaxa = hiddenState.getRootedGeneGenealogy().getLeaves();
		Arrays.sort(sortedCheckTaxa);
		if (!Arrays.equals(sortedReferenceTaxa, sortedCheckTaxa)) {
		    System.err.println ("ERROR: list of taxa in allele-species-mapping is not consistent with gene genealogies.");
		    System.err.println (Arrays.toString(sortedReferenceTaxa) + "\n");
		    System.err.println (Arrays.toString(sortedCheckTaxa) + "\n");
		    return (false);
		}
	    }
	}

	return (true);
    }

    // bleh - special case for outgroup taxon
    // always mapped to same name for species
    // ignore this
    protected boolean verifyParentalTreeTaxa () {
	// check parental tree taxa vs. species-alleles-mapping keys
	for (HiddenState hiddenState : hiddenStates) {
	    String[] sortedReferenceTaxa = getSortedTaxaFromSpeciesToAllelesMapping(hiddenState.getSpeciesToAllelesMapping(), false);
	    String[] sortedCheckTaxa = getTaxa(hiddenState.getParentalTree());
	    Arrays.sort(sortedCheckTaxa);
	    if (!Arrays.equals(sortedReferenceTaxa, sortedCheckTaxa)) {
		System.err.println ("ERROR: list of taxa in basic-info-file is not consistent with parental trees.");
		System.err.println (Arrays.toString(sortedReferenceTaxa) + "\n");
		System.err.println (Arrays.toString(sortedCheckTaxa) + "\n");
		return (false);
	    }
	}

	return (true);
    }

    /** 
     * Set geneGenealogyTaxaFlag to true to check rooted gene genealogy taxa, set to false
     * to check parental tree taxa.
     *
     * Unrooted gene genealogy taxa match rooted gene genealogy taxa by construction.
     */
    // protected boolean verifyTaxa (String[] sortedReferenceTaxa, boolean geneGenealogyTaxaFlag) {
    // 	for (HiddenState hiddenState : hiddenStates) {
    // 	    String[] sortedCheckTaxa = geneGenealogyTaxaFlag ? 
    // 		hiddenState.getRootedGeneGenealogy().getLeaves() : 
    // 		getTaxa(hiddenState.getParentalTree());
    // 	    Arrays.sort(sortedCheckTaxa);
    // 	    if (!Arrays.equals(sortedReferenceTaxa, sortedCheckTaxa)) {
    // 		return (false);
    // 	    }
    // 	}

    // 	return (true);
    // }

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
    protected String[] getTaxa (Network<CoalescePattern[]> network) {
	Vector<String> names = new Vector<String>();
	for (NetNode<CoalescePattern[]> leaf : network.getLeaves()) {
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
	    //System.out.println("2) Learn Model using Grid Search");
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
    // protected int sequenceChoice(BufferedReader in) {
    // 	boolean init = true;
    // 	int option = -1;
    // 	while (init) {
    // 	    System.out.println("Which observation sequence would you like to use?");
    // 	    System.out.println("0) Reuse previously read sequence.");
    // 	    System.out.println("1) Load a new observation sequence.");
    // 	    System.out.println("Choose an option: ");
    // 	    option = getOption(2, in);
    // 	    if (option != -1) init = false;
    // 	}
    // 	return option;
    // }


    /**
     * Read and Get Observation
     */
    protected ArrayList<ArrayList<ObservationMap>> getDataPartitions(BufferedReader in) {
	try {
	    System.out.println("Keep or discard parsimony-uninformative sites? true for keep, false for discard: ");
	    boolean keepUninformativeSitesFlag = Boolean.parseBoolean(in.readLine());
	    System.out.println("Input whitespace-delimited list of FASTA alignments: ");
	    String filenamesString = in.readLine();
	    StringTokenizer st = new StringTokenizer(filenamesString);
	    ArrayList<ArrayList<ObservationMap>> dataPartitions = new ArrayList<ArrayList<ObservationMap>>();
	    while (st.hasMoreTokens()) {
		String filename = st.nextToken();
		fParser.parseMe(filename, keepUninformativeSitesFlag);
		dataPartitions.add(fParser.getObs());
	    }
	    return (dataPartitions);
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
     * 
     * No switching frequency ratio term contribution.
     */
    protected double[] calculatePi () {
	double[] pi = new double[hiddenStates.size()];
	double norm = 0.0;
	for (int i = 0; i < hiddenStates.size(); i++) {
	    HiddenState hiddenState = hiddenStates.get(i);
	    pi[i] = hiddenState.calculateProbabilityOfRootedGeneGenealogyInParentalTree();
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

    /**
     * Do two hidden states share the same parental tree object?
     * (Same row in the HMM hidden state space?)
     */
    protected boolean checkSameParentalTreeClass (HiddenState si, HiddenState sj) {
	Set<HiddenState> sic = parentalTreeClasses.get(si.getParentalTree());
	return (sic.contains(sj));
    }

    /**
     * Do two hidden states share the same gene genealogy object?
     * (Same column in the HMM hidden state space?)
     */
    protected boolean checkSameGeneRootedGenealogyClass (HiddenState si, HiddenState sj) {
	Set<HiddenState> sic = rootedGeneGenealogyClasses.get(si.getRootedGeneGenealogy());
	return (sic.contains(sj));
    }

    protected boolean checkSameGeneUnrootedGenealogyClass (HiddenState si, HiddenState sj) {
	Set<HiddenState> sic = unrootedGeneGenealogyClasses.get(si.getUnrootedGeneGenealogy());
	return (sic.contains(sj));
    }

    // Set<HiddenState> namesSet, BidirectionalMultimap<Tuple<HiddenState,HiddenState>,SwitchingFrequencyRatioTerm> pairToSwitchingFrequencyRatioTermMap, 

    /**
     * Solve for switching frequencies based on switching frequency ratio terms.
     * Set parentalTreeFlag to true to calculate for parental-tree-switching,
     * or false to calculate for gene-genealogy-switching.
     *
     * Don't constrain switching frequencies to be at most 1/k, for k the number of alternatives.
     * Might be too small of a value.
     * Instead, if self-switching frequency falls below 1/(1 + SwitchingFrequencyRatioTermParameter.DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_TERMS_TOTAL_WEIGHT), 
     * then normalize+rescale non-self-transition-switching-frequency-ratio-terms 
     * to have total weight of DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_TERMS_TOTAL_WEIGHT
     * (leaving self-transition-frequency-ratio-term as 1).
     * Kind of weird non-linear behavior, but oh well. We expect funky optimization landscape anyways.
     */
    // protected MapOfMap<HiddenState,HiddenState,Double> calculateSwitchingFrequencies (boolean verifyFlag) {
    // 	MapOfMap<HiddenState,HiddenState,Double> result = new MapOfMap<HiddenState,HiddenState,Double>();
    // 	for (int i = 0; i < hiddenStates.size(); i++) {
    // 	    double totalWeight = 0.0;
    // 	    for (int j = 0; j < hiddenStates.size(); j++) {
    // 		if (i == j) {
    // 		    // canonically, self-switching frequency has ratio term 1.0
    // 		    result.put(hiddenStates.get(i), hiddenStates.get(j), new Double(1.0));
    // 		    totalWeight += 1.0;
    // 		    continue;
    // 		}

    // 		// bleh - final
    // 		Tuple<HiddenState,HiddenState> pair = new Tuple<HiddenState,HiddenState>(hiddenStates.get(i), hiddenStates.get(j));
    // 		Set<SwitchingFrequencyRatioTerm> terms = hiddenStatePairToSwitchingFrequencyRatioTermMap.get(pair);
    // 		// paranoid
    // 		if ((terms == null) || (terms.size() <= 0) || (terms.size() > 1)) {
    // 		    System.err.println("ERROR: unable to retrieve switching-frequency-ratio-term for named pair " + pair + ".");
    // 		    return (null);
    // 		}
		
    // 		SwitchingFrequencyRatioTerm term = terms.iterator().next();

    // 		// normalize by total weight
    // 		result.put(hiddenStates.get(i), hiddenStates.get(j), new Double(term.getValue()));
    // 		totalWeight += term.getValue();
    // 	    }

    // 	    // paranoid
    // 	    if (totalWeight <= 0.0) {
    // 		throw (new RuntimeException("ERROR: totalWeight is non-positive in runHmm.calculateSwitchingFrequencies(). " + hiddenStates.get(i).getName() + " " + totalWeight));
    // 	    }

    // 	    // no - causes seriously broken behavior
    // 	    // forces any self-transition to have minimum 2/3 probability
    // 	    // even for a gene genealogy that is incongruent within its containing parental tree

    // 	    // // if total weight on non-self-transitions is too much, 
    // 	    // // normalize+rescale non-self-transitions' switching-frequency-ratio-terms
    // 	    // if (totalWeight - 1.0 >= optimize.SwitchingFrequencyRatioTermParameter.DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_TERMS_TOTAL_WEIGHT) {
    // 	    // 	if (Constants.WARNLEVEL > 4) { System.out.println ("INFO: total weight of non-self-transitions away from hidden state " + hiddenStates.get(i).getName() + " is larger than " + optimize.SwitchingFrequencyRatioTermParameter.DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_TERMS_TOTAL_WEIGHT + ". Normalizing and re-scaling non-self-transitions for this hidden state."); }

    // 	    // 	double updatedTotalWeight = 1.0;
    // 	    // 	for (int j = 0; j < hiddenStates.size(); j++) {
    // 	    // 	    if (i == j) {
    // 	    // 		continue;
    // 	    // 	    }

    // 	    // 	    // wonky - doesn't update the SwitchingFrequencyRatioTerm itself
    // 	    // 	    // bizarre behavior - same parameter can "saturate" for one row but not another
    // 	    // 	    double originalWeight = result.get(hiddenStates.get(i), hiddenStates.get(j)).doubleValue();
    // 	    // 	    double newWeight = originalWeight / ((totalWeight - 1.0) / optimize.SwitchingFrequencyRatioTermParameter.DEFAULT_MAXIMUM_TOTAL_NON_SELF_TRANSITION_TERMS_TOTAL_WEIGHT);
    // 	    // 	    result.put(hiddenStates.get(i), hiddenStates.get(j), new Double(newWeight));
    // 	    // 	    updatedTotalWeight += newWeight;
    // 	    // 	}
    // 	    // 	totalWeight = updatedTotalWeight;
    // 	    // }

    // 	    // now normalize all transitions away from hidden state i
    // 	    // to get probabilities
    // 	    for (int j = 0; j < hiddenStates.size(); j++) {
    // 		double originalWeight = result.get(hiddenStates.get(i), hiddenStates.get(j)).doubleValue();
    // 		result.put(hiddenStates.get(i), hiddenStates.get(j), new Double(originalWeight / totalWeight));
    // 	    }
    // 	}

    // 	// paranoid
    // 	if (verifyFlag) {
    // 	    for (int i = 0; i < hiddenStates.size(); i++) {
    // 		double totalWeight = 0.0;
    // 		//double[] row = new double[hiddenStates.size()];
    // 		for (int j = 0; j < hiddenStates.size(); j++) {
    // 		    double fij = result.get(hiddenStates.get(i), hiddenStates.get(j));
    // 		    if ((fij < 0.0) || (fij > 1.0)) {
    // 			System.err.println ("ERROR: switching-frequency for named pair " + hiddenStates.get(i).getName() + " " + hiddenStates.get(j).getName() + " isn't a probability. Returning null to signal error.");
    // 			return (null);
    // 		    }
    // 		    totalWeight += fij;
    // 		    //row[j] = fij;
    // 		}
    // 		if (Math.abs(totalWeight - 1.0) > Constants.ZERO_DELTA) {
    // 		    //System.out.println ("offending row: |" + Arrays.toString(row) + "|");

    // 		    System.err.println ("ERROR: row-summed-switching-frequencies for row named " + hiddenStates.get(i).getName() + " doesn't sum to one. Returning null to signal error. " + totalWeight);
    // 		    return (null);
    // 		}
    // 	    }
    // 	}

    // 	return (result);
    // }

    // boolean verifyFlag
    
    // one-shot it
    // function should be renamed - really calculating full transition probabilities in one shot
    protected MapOfMap<HiddenState,HiddenState,Double> calculateSwitchingFrequencies () {
	MapOfMap<HiddenState,HiddenState,Double> result = new MapOfMap<HiddenState,HiddenState,Double>();

	// calculate within-row transition probabilities
	for (String p1name : parentalTreeNameMap.keys()) {
	    for (String g1name : geneGenealogyNameMap.keys()) {
		HiddenState h1 = parentalTreeGeneGenealogyNamePairToHiddenStateMap.get(p1name, g1name);
		double totalWeight = 0.0;
		for (String g2name : geneGenealogyNameMap.keys()) {
		    HiddenState h2 = parentalTreeGeneGenealogyNamePairToHiddenStateMap.get(p1name, g2name);
		    // if (g1name.equals(g2name)) {
		    // 	result.put(h1, h2, new Double(1.0));
		    // 	totalWeight += 1.0;
		    // }
		    // else {
		    Tuple<HiddenState,HiddenState> pair = new Tuple<HiddenState,HiddenState>(h1, h2);
		    Set<SwitchingFrequencyRatioTerm> terms = hiddenStatePairToSwitchingFrequencyRatioTermMap.get(pair);
		    // paranoid
		    if ((terms == null) || (terms.size() <= 0) || (terms.size() > 1)) {
			System.err.println("ERROR: unable to retrieve switching-frequency-ratio-term for named pair " + pair + ".");
			return (null);
		    }
		
		    SwitchingFrequencyRatioTerm term = terms.iterator().next();
		    result.put(h1, h2, new Double(term.getValue()));
		    totalWeight += term.getValue();
			//}
		}

		// convert ratios to frequencies 
		for (String g2name : geneGenealogyNameMap.keys()) {
		    HiddenState h2 = parentalTreeGeneGenealogyNamePairToHiddenStateMap.get(p1name, g2name);
		    double originalWeight = result.get(h1, h2).doubleValue();
		    result.put(h1, h2, new Double(originalWeight / totalWeight));
		}

		double totalWeight2 = 0.0;
		for (String g2name : geneGenealogyNameMap.keys()) {
		    HiddenState h2 = parentalTreeGeneGenealogyNamePairToHiddenStateMap.get(p1name, g2name);
		    result.put (h1, h2, new Double(
						   result.get(h1, h2).doubleValue() * 
						   h2.calculateProbabilityOfRootedGeneGenealogyInParentalTree()
						   )
				);
		    // order is critical
		    totalWeight2 += result.get(h1, h2).doubleValue();
		}

		// now normalize within-row, and also scale by (1-\gamma)
		for (String g2name : geneGenealogyNameMap.keys()) {
		    HiddenState h2 = parentalTreeGeneGenealogyNamePairToHiddenStateMap.get(p1name, g2name);
		    double originalWeight = result.get(h1, h2).doubleValue();
		    result.put(h1, h2, new Double((originalWeight / totalWeight2) * (1.0 - gamma.getValue())));
		}
	    }
	}

	// now calculate across-row switching frequencies
	for (String p1name : parentalTreeNameMap.keys()) {
	    for (String p2name : parentalTreeNameMap.keys()) {
		if (p1name.equals(p2name)) {
		    continue;
		}

		for (String g1name : geneGenealogyNameMap.keys()) {
		    HiddenState h1 = parentalTreeGeneGenealogyNamePairToHiddenStateMap.get(p1name, g1name);
		    for (String g2name : geneGenealogyNameMap.keys()) {
			HiddenState h2 = parentalTreeGeneGenealogyNamePairToHiddenStateMap.get(p2name, g2name);

			// only \gamma and P[g|T] ILS contribution for across-row switching
			// keep it simple for now
			//
			// if this looks promising, can try adding more 
			// parameterization here, as needed
			result.put(h1, h2, new Double(
						      h2.calculateProbabilityOfRootedGeneGenealogyInParentalTree() *
						      gamma.getValue() / ((gamma.getNumAlternatives() - 1) * 1.0)
						      ));
		    }
		}
	    }
	}

	return (result);
    }


    /**
     * Calculate initial transition probability matrix a_{ij}.
     * See revised writeup for details.
     *
     * Call this after changing parental tree branch lengths, hybridization frequency, and
     * any other parameters related to the transition probabilities.
     */
    protected double[][] calculateAij () {
	// may want to cache 
	// parental-tree-switching and gene-genealogy-switching probabilities
	//
	// Don't get too fancy. Just base switching on names of parental trees and gene genealogies.
	//
	// parental-tree-switching frequencies, computed using parental-tree-switching frequency ratio terms
	//
	// if no cache entry, re-calculate
	// if ((calculationCache.cacheSwitchingFrequencyMap == null) || calculationCache.cacheSwitchingFrequencyMap.isEmpty()) {
	//     MapOfMap<HiddenState,HiddenState,Double> switchingFrequencyMap = calculateSwitchingFrequencies(true);
	//     calculationCache.cacheSwitchingFrequencyMap = switchingFrequencyMap;
	// }
	// MapOfMap<HiddenState,HiddenState,Double> switchingFrequencyMap = calculationCache.cacheSwitchingFrequencyMap;

	MapOfMap<HiddenState,HiddenState,Double> switchingFrequencyMap = calculateSwitchingFrequencies();

	double[][] a = new double[hiddenStates.size()][hiddenStates.size()];
	for (int i = 0; i < a.length; i++) {
	    HiddenState si = hiddenStates.get(i);
	    //String pi = parentalTreeNameMap.rget(si.getParentalTree());
	    //String gi = geneGenealogyNameMap.rget(si.getRootedGeneGenealogy());
	    for (int j = 0 ; j < a[i].length; j++) {
		HiddenState sj = hiddenStates.get(j);
		//String pj = parentalTreeNameMap.rget(sj.getParentalTree());
		//String gj = geneGenealogyNameMap.rget(sj.getRootedGeneGenealogy());
		a[i][j] = 
		    switchingFrequencyMap.get(si, sj).doubleValue()
		    ;
	    }
	}

	// bleh - post new model extensions,
	// need to re-normalize after (parental-tree-switching parameter) x (gene-genealogy-switching parameter) x
	// (ILS calculation from Degnan and Salter 2005) contributions all factored in
	//rowNormalize(a);

	// strict!
	if (!verifyAij(a)) {
	    System.err.println ("ERROR: verifyAij() failed. Returning null to signal error.");
	    return (null);
	}

	return (a);
    }    

    protected void rowNormalize (double[][] a) {
	double norm;
	for (int i = 0; i < a.length; i++) {
	    norm = 0.0;
	    for (int j = 0; j < a[i].length; j++) {
		norm += a[i][j];
	    }

	    // paranoid
	    if (norm <= 0.0) {
		throw (new RuntimeException ("ERROR: norm is non-positive in rowNormalize()."));
	    }

	    for (int j = 0; j < a[i].length; j++) {
		a[i][j] /= norm;
	    }
	}
    }

    /**
     * Verify that each entry is a probability. 
     * Also verify that each row sums to one.
     */
    protected boolean verifyAij (double[][] a) {	
	for (int i = 0; i < a.length; i++) {
	    double rowsum = 0.0;
	    for (int j = 0; j < a[i].length; j++) {
		if ((a[i][j] < 0.0) || (a[i][j] > 1.0)) {
		    System.err.println ("ERROR: entry in a_ij transition matrix is invalid. " + a[i][j]);
		    return (false);
		}
		rowsum += a[i][j];

		//if (Constants.WARNLEVEL > 4) { System.out.println(a[i][j] + " "); }
	    }
	    
	    //if (Constants.WARNLEVEL > 4) { System.out.println(); }

	    if (Math.abs(rowsum - 1.0) > Constants.ZERO_DELTA) {
		System.err.println ("ERROR: row " + i + " in a_ij transition matrix doesn't sum to one. Sum is " + rowsum + ". Make sure that list of gene genealogies includes all possible topologies.");
		return (false);
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
	// kliu - gene genealogy branch lengths are getting lost in here
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

	myhmm = buildMyHmm(hiddenStates, pi, a);
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
	fParser.setTrees(hiddenStates);	
    }
	
    

    // EvoTree etree
    // etree.toNewickString(true, true)
    /**
     * Convert an EvoTree into a PhyloNet Network object
     */
    protected Network<CoalescePattern[]> convertPHMMTreeToPhyloNetNetwork (String newickString) {
    	ExNewickReader<CoalescePattern[]> enr = new ExNewickReader<CoalescePattern[]>(new StringReader(newickString));
    	// kliu - change this over to a String data member
    	// don't really use it anyways
    	try {
    	    Network<CoalescePattern[]> network = enr.readNetwork();
    	    return (network);
    	}
    	catch(IOException ioe) {
    	    System.err.println(ioe);
    	    ioe.printStackTrace();
    	    return (null);
    	}
    }

    // EvoTree etree
    // etree.toNewickString(true, true)
    /**
     * Convert an EvoTree into a PhyloNet Tree object.
     */
    protected Tree convertNewickStringToPhyloNetTree (String newickString) {
	// kliu - weird - this is losing the branch length information?
    	NewickReader nr = new NewickReader(new StringReader(newickString));
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
     * Expose this so other classes can access the class name for an unrooted gene genealogy.
     */
    public String getTopologicalEquivalenceClassName (Tree unrootedGeneGenealogy) {
	Vector<String> names = new Vector<String>();
	for (Tree rootedMember : rootedToUnrootedGeneGenealogyMap.rget(unrootedGeneGenealogy)) {
	    names.add(geneGenealogyNameMap.rget(rootedMember));
	}

	Collections.sort(names);

	String result = "";
	for (int i = 0; i < names.size(); i++) {
	    if (i > 0) {
		result += HiddenState.EQUIVALENCE_CLASS_NAME_DELIMITER;
	    }
	    result += names.get(i);
	}

	return (result);
    }

    /**
     * Read and parse trees file
     * @throws Exception 
     */
    protected void buildTrees() throws Exception {
	if ((parentalTreesFileName == null) ||
	    (!(new File(parentalTreesFileName)).exists())) {
	    throw new ParserFileException("ERROR: invalid parental trees file.");
	}

	hiddenStates = new ArrayList<HiddenState>();
	// also maintain equivalence classes among hidden states based on shared parental trees
	parentalTreeClasses = new Hashtable<Network<CoalescePattern[]>,Set<HiddenState>>();
	rootedGeneGenealogyClasses = new Hashtable<Tree,Set<HiddenState>>();
	unrootedGeneGenealogyClasses = new Hashtable<Tree,Set<HiddenState>>();
	// also retain parental tree names
	// to facilitate parameter inputs/constraints on parental tree branches
	parentalTreeNameMap = new BijectiveHashtable<String,Network<CoalescePattern[]>>();
	geneGenealogyNameMap = new BijectiveHashtable<String,Tree>();
	// retain for construction of HiddenState objects
	parentalTreeSpeciesToAllelesMapMap = new Hashtable<String,Map<String,List<String>>>();
	parentalTreeGeneGenealogyNamePairToHiddenStateMap = new MapOfMap<String,String,HiddenState>();

	System.out.println("\nNow building trees . . .");
	BufferedReader ptreesbr = new BufferedReader(new FileReader(parentalTreesFileName));
	String line;
	while ((line = ptreesbr.readLine()) != null) {
	    StringTokenizer st = new StringTokenizer(line);
	    if (st.countTokens() != 3) {
		throw (new RuntimeException("ERROR: invalid parental trees line " + line + " in file " + parentalTreesFileName + "."));
	    }
	    String parentalTreeName = st.nextToken();
	    Network<CoalescePattern[]> parentalTree = convertPHMMTreeToPhyloNetNetwork(st.nextToken());
	    Map<String,List<String>> speciesToAllelesMap = parseSpeciesToAllelesMap(st.nextToken());
	    // no duplicate node names allowed!
	    if (parentalTree.hasDuplicateNames()) {
		throw (new RuntimeException("ERROR: duplicate node names are present in parental tree " + parentalTreeName + ". Check inputs and try again."));
	    }
	    // no duplicate parental tree names allowed!
	    if (parentalTreeNameMap.containsKey(parentalTreeName)) {
		throw (new RuntimeException("ERROR: duplicate parental tree name " + parentalTreeName + ". Check inputs and try again."));
	    }
	    parentalTreeNameMap.put(parentalTreeName, parentalTree);
	    parentalTreeSpeciesToAllelesMapMap.put(parentalTreeName, speciesToAllelesMap);
	}
	ptreesbr.close();
	
	BufferedReader ggbr = new BufferedReader(new FileReader(geneGenealogiesFileName));
	while ((line = ggbr.readLine()) != null) {
	    StringTokenizer st = new StringTokenizer(line);
	    if (st.countTokens() != 2) {
		throw (new RuntimeException("ERROR: invalid gene genealogies line " + line + " in file " + geneGenealogiesFileName + "."));
	    }
	    String geneGenealogyName = st.nextToken();
	    Tree geneGenealogy = convertNewickStringToPhyloNetTree(st.nextToken());
	    if (hasDuplicateNames(geneGenealogy)) {
		throw (new RuntimeException("ERROR: duplicate node names are present in gene genealogy " + geneGenealogyName + ".Check inputs and try again."));
	    }

	    // really strict
	    if (parentalTreeNameMap.containsKey(geneGenealogyName)) {
		throw (new RuntimeException("ERROR: gene genealogy name " + geneGenealogyName + " appears as a parental tree name. Check inputs and try again."));
	    }
	    
	    if (geneGenealogyNameMap.containsKey(geneGenealogyName)) {
		throw (new RuntimeException("ERROR: duplicate gene genealogy name " + geneGenealogyName + ". Check inputs and try again."));
	    }
	    
	    geneGenealogyNameMap.put(geneGenealogyName, geneGenealogy);	    
	}
	ggbr.close();

	// Now create equivalence classes among rooted gene genealogies based on 
	// unrooted topological equality.
	//
	// Only optimize branch lengths of unrooted member from each topological equivalence
	// class. To prevent overfitting.
	rootedToUnrootedGeneGenealogyMap = createRootedToUnrootedGeneGenealogyMap();

	// testing
	// if (Constants.WARNLEVEL > 2) { System.out.println ("Collapse gene genealogies by topological equivalence class? " + collapseGeneGenealogiesByTopologicalEquivalenceClassFlag); System.out.println(); }
	if (Constants.WARNLEVEL > 2) { debugRootedToUnrootedGeneGenealogyMap(); }

	// kliu - indexing is by (parentalTree, geneGenealogy) appearance order according to the following:
	// kliu - hmm... would be nice to go in alphabetical order according to tree names
	ArrayList<String> parentalTreeNames = new ArrayList<String>(parentalTreeNameMap.keys());
	Collections.sort(parentalTreeNames);
	ArrayList<String> geneGenealogyNames = new ArrayList<String>(geneGenealogyNameMap.keys());
	Collections.sort(geneGenealogyNames);
	for (String parentalTreeName : parentalTreeNames) {
	    Network<CoalescePattern[]> parentalTree = parentalTreeNameMap.get(parentalTreeName);
	    for (String geneGenealogyName : geneGenealogyNames) {
		Tree rootedGeneGenealogy = geneGenealogyNameMap.get(geneGenealogyName);
		// paranoid
		if (!rootedToUnrootedGeneGenealogyMap.containsKey(rootedGeneGenealogy) || 
		    rootedToUnrootedGeneGenealogyMap.get(rootedGeneGenealogy).size() != 1) {
		    throw (new RuntimeException("ERROR: invalid entry for rooted gene genealogy " + geneGenealogyName + " in rootedToUnrootedGeneGenealogyMap in runHmm.buildTrees(...)."));
		}
		Tree unrootedGeneGenealogy = rootedToUnrootedGeneGenealogyMap.get(rootedGeneGenealogy).iterator().next();
		String hiddenStateName = parentalTreeNameMap.rget(parentalTree) + HiddenState.HIDDEN_STATE_NAME_DELIMITER + geneGenealogyNameMap.rget(rootedGeneGenealogy);
		//if (collapseGeneGenealogiesByTopologicalEquivalenceClassFlag) {
		    hiddenStateName += HiddenState.HIDDEN_STATE_NAME_DELIMITER + getTopologicalEquivalenceClassName(unrootedGeneGenealogy);
		    //}
		// kliu - meh - parse allele-to-species mapping later and add in references here
		HiddenState hiddenState = new HiddenState(hiddenStateName, 
							  parentalTree, 
							  rootedGeneGenealogy, 
							  //collapseGeneGenealogiesByTopologicalEquivalenceClassFlag ? unrootedGeneGenealogy : rootedGeneGenealogy, 
							  unrootedGeneGenealogy,
							  outgroupTaxonName,
							  parentalTreeSpeciesToAllelesMapMap.get(parentalTreeName), 
							  gtrSubstitutionModel, 
							  calculationCache);
		hiddenStates.add(hiddenState);
		// uber paranoid - unneeded due to name uniqueness of parental trees and gene genealogies,
		// but oh well
		if (parentalTreeGeneGenealogyNamePairToHiddenStateMap.contains(parentalTreeName, geneGenealogyName)) {
		    throw (new RuntimeException("ERROR: hidden state already exists for tree pair named " + parentalTreeName + " " + geneGenealogyName + "."));
		}
		parentalTreeGeneGenealogyNamePairToHiddenStateMap.put(parentalTreeName, geneGenealogyName, hiddenState);
		

		// maintain equivalence class maps
		if (!parentalTreeClasses.containsKey(parentalTree)) {
		    parentalTreeClasses.put(parentalTree, new HashSet<HiddenState>());
		}
		
		if (!rootedGeneGenealogyClasses.containsKey(rootedGeneGenealogy)) {
		    rootedGeneGenealogyClasses.put(rootedGeneGenealogy, new HashSet<HiddenState>());
		}

		if (!unrootedGeneGenealogyClasses.containsKey(unrootedGeneGenealogy)) {
		    unrootedGeneGenealogyClasses.put(unrootedGeneGenealogy, new HashSet<HiddenState>());
		}

		parentalTreeClasses.get(parentalTree).add(hiddenState);
		rootedGeneGenealogyClasses.get(rootedGeneGenealogy).add(hiddenState);
		unrootedGeneGenealogyClasses.get(unrootedGeneGenealogy).add(hiddenState);
	    }
	}
	
	// ------- >Testing purposes
	//	        System.out.println("Trees read in and built:");
	//	        for (int i = 0; i < hiddenStates.size(); i++) {
	//	        	System.out.println((hiddenStates.get(i)));
	//	        }
	// <------Testing purposes

	if (numStates != hiddenStates.size()) {
	    throw new ParserFileException("Error: the number of trees read is not the same as the number of states previously inputted!");
	}	

	// strict!
	if (hiddenStates.size() < 2) {
	    throw (new ParserFileException("ERROR: must be at least two hidden states in PhyloNet-HMM."));
	}
    }

    /**
     * kliu - reverse the mapping per the latest PhyloNet changes.
     */
    protected Map<String,List<String>> parseSpeciesToAllelesMap (String mapString) {
	Map<String,List<String>> map = new Hashtable<String,List<String>>();
	StringTokenizer st1 = new StringTokenizer(mapString, SPECIES_TO_ALLELES_MAPPING_ENTRY_DELIMITER);
	while (st1.hasMoreTokens()) {
	    StringTokenizer st2 = new StringTokenizer(st1.nextToken(), SPECIES_TO_ALLELES_MAPPING_KEY_VALUE_DELIMITER);
	    if (st2.countTokens() != 2) {
		// strict!
		throw (new RuntimeException("ERROR: incorrectly formatted allele-to-species map. " + mapString));
	    }

	    String allele = st2.nextToken();
	    String species = st2.nextToken();
	    if (!map.containsKey(species)) {
		map.put(species, new Vector<String>());
	    }

	    map.get(species).add(allele);
	}
	return (map);
    }

    /**
     * Print out equivalence classes.
     */
    protected void debugRootedToUnrootedGeneGenealogyMap () {
	System.out.println ("Gene genealogy topological equivalence classes: ");
	for (Tree unrootedGeneGenealogy : rootedToUnrootedGeneGenealogyMap.values()) {
	    System.out.println ("Unrooted gene genealogy " + unrootedGeneGenealogy.toNewick() + " topological equivalence class members: ");
	    for (Tree rootedGeneGenealogy : rootedToUnrootedGeneGenealogyMap.rget(unrootedGeneGenealogy)) {
		System.out.println (rootedGeneGenealogy.toNewick());
	    }
	    System.out.println();
	}
	System.out.println();
    }

    /**
     * Inefficient O(n^3) code. Bleh.
     * But this function only gets run once during an analysis.
     * Don't spend time being clever about this function.
     *
     * We optimize only the unrooted representative (whose rooting has no meaning) of a topological equivalence class.
     */
    protected BidirectionalMultimap<Tree,Tree> createRootedToUnrootedGeneGenealogyMap () {
	// paranoid
	if ((geneGenealogyNameMap == null) || (geneGenealogyNameMap.sizeKeys() <= 0)) {
	    System.err.println ("ERROR: called runHmm.createRootedToUnrootedGeneGenealogyMap() with null or empty geneGenealogyNameMap. Returning null to signal error.");
	    return (null);
	}

	Tree[] geneGenealogies = new Tree[geneGenealogyNameMap.values().size()];
	geneGenealogies = geneGenealogyNameMap.values().toArray(geneGenealogies);

	// bleh
	// geneGenealogies index -> equivalence class
	// weird Vector behavior
	// need to call setSize() to force Vector to take on a fixed size - bleh
	Vector<Set<Tree>> equivalenceClassMap = new Vector<Set<Tree>>();
	// weird Vector behavior with initial capacity
	// ArrayList doesn't even support setSize()
	equivalenceClassMap.setSize(geneGenealogies.length); 
	

	// not the most efficient way to compute pairwise R-F distance matrix
	// Day's algorithm gives a linear time algo for a set of trees
	// meh, still need to read/write n^2 trees just for matrix itself
	for (int i = 0; i < geneGenealogies.length; i++) {
	    for (int j = i + 1; j < geneGenealogies.length; j++) {
		// kliu - seems to be failing for trees that lack internal edges
		// same class
		if (TreeUtils.calculateRobinsonFouldsDistance(geneGenealogies[i], geneGenealogies[j]) == 0) {
		    if ((equivalenceClassMap.get(i) != null) &&
			(equivalenceClassMap.get(j) == null)) {
			equivalenceClassMap.get(i).add(geneGenealogies[j]);
			equivalenceClassMap.set(j, equivalenceClassMap.get(i));
		    }
		    else if ((equivalenceClassMap.get(i) == null) &&
			     (equivalenceClassMap.get(j) != null)) {
			equivalenceClassMap.get(j).add(geneGenealogies[i]);
			equivalenceClassMap.set(i, equivalenceClassMap.get(j));
		    }
		    else if ((equivalenceClassMap.get(i) == null) &&
			     (equivalenceClassMap.get(j) == null)) {
			Set<Tree> st = new HashSet<Tree>();
			st.add(geneGenealogies[i]);
			st.add(geneGenealogies[j]);
			equivalenceClassMap.set(i, st);
			equivalenceClassMap.set(j, st);
		    }
		    else { // both sets exist
			// coalesce
			// expensive op
			// makes this function O(n^3)
			// bleh
			// later, spend time to make this n^2?
			Set<Tree> smallerSet = equivalenceClassMap.get(i);
			Set<Tree> largerSet = equivalenceClassMap.get(j);
			if (equivalenceClassMap.get(j).size() < equivalenceClassMap.get(i).size()) {
			    smallerSet = equivalenceClassMap.get(j);
			    largerSet = equivalenceClassMap.get(i);
			}
			largerSet.addAll(smallerSet);
			// rather than keep track of i vs. j
			// just do a repeated set
			// it's fine
			equivalenceClassMap.set(i, largerSet);
			equivalenceClassMap.set(j, largerSet);
		    }
		}
		else { // different class
		    if (equivalenceClassMap.get(i) == null) {
			Set<Tree> st = new HashSet<Tree>();
			st.add(geneGenealogies[i]);
			equivalenceClassMap.set(i, st);
		    }

		    if (equivalenceClassMap.get(j) == null) {
			Set<Tree> st = new HashSet<Tree>();
			st.add(geneGenealogies[j]);
			equivalenceClassMap.set(j, st);
		    }
		}
	    }
	}

	// now get unique classes from map
	Set<Set<Tree>> uniqueEquivalenceClasses = new HashSet<Set<Tree>>();
	for (Set<Tree> equivalenceClass : equivalenceClassMap) {
	    if (!uniqueEquivalenceClasses.contains(equivalenceClass)) {
		uniqueEquivalenceClasses.add(equivalenceClass);
	    }
	}

	BidirectionalMultimap<Tree,Tree> result = new BidirectionalMultimap<Tree,Tree>();

	for (Set<Tree> equivalenceClass : uniqueEquivalenceClasses) {
	    Tree canonicalRepresentative = computeCanonicalRepresentativeOfTopologicalEquivalenceClass(equivalenceClass);
	    // create a new copy
	    // looks a little nicer if input rooted gene genealogies are left intact
	    Tree canonicalRepresentativeCopy = new STITree<Double>(canonicalRepresentative);
	    for (Tree member : equivalenceClass) {
		result.put(member, canonicalRepresentativeCopy);
	    }
	}

	return (result);
    }

    /**
     * The rooting of the canonical representative is totally irrelevant,
     * since we're using reversible substitution models.
     * The branch lengths shouldn't make a difference if the optimization heuristic
     * is decent.
     * Only do it this way since Felsenstein code uses a rooted tree data structure.
     *
     * Just want to optimize a single representative of a topological equivalence class,
     * to prevent over-fitting.
     */
    protected Tree computeCanonicalRepresentativeOfTopologicalEquivalenceClass (Set<Tree> equivalenceClass) {
	if (equivalenceClass.size() <= 0) {
	    throw (new RuntimeException ("ERROR: empty equivalence class in runHmm.computeCanonicalRepresentativeOfTopologicalEquivalenceClass(...)."));
	}
	return (equivalenceClass.iterator().next());
    }

    // parentalTreeEquivalenceClass
	// TreeParser ptp = new TreeParser(ptreesbr);
	// ArrayList<EvoTree> eParentalTrees = ptp.nexusFileTreeNames(parentalTreesFileName);
	// for (EvoTree etree : eParentalTrees) {
	//     Network<CoalescePattern[]> parentalTree = convertPHMMTreeToPhyloNetNetwork(etree);
	//     // no duplicate node names allowed!
	//     if (parentalTree.hasDuplicateNames()) {
	// 	throw (new RuntimeException("ERROR: duplicate node names are present in parental tree " + etree.getName() + ". Check inputs and try again."));
	//     }
	//     // no duplicate parental tree names allowed!
	//     if (parentalTreeNameMap.containsKey(etree.getName())) {
	// 	throw (new RuntimeException("ERROR: duplicate parental tree name " + etree.getName() + ". Check inputs and try again."));
	//     }
	//     parentalTreeNameMap.put(etree.getName(), parentalTree);
	// }
	// ptreesbr.close();
    
	// TreeParser gtp = new TreeParser(ggbr);
	// ArrayList<EvoTree> eGeneGenealogies = gtp.nexusFileTreeNames(geneGenealogiesFileName);
	// ggbr.close();

	// for (EvoTree egg : eGeneGenealogies) {
	//     Tree geneGenealogy = convertPHMMTreeToPhyloNetTree(egg);
	//     if (hasDuplicateNames(geneGenealogy)) {
	// 	throw (new RuntimeException("ERROR: duplicate node names are present in gene genealogy " + egg.getName() + ".Check inputs and try again."));
	//     }

	//     // really strict
	//     if (parentalTreeNameMap.containsKey(egg.getName())) {
	// 	throw (new RuntimeException("ERROR: gene genealogy name " + egg.getName() + " appears as a parental tree name. Check inputs and try again."));
	//     }
	    
	//     if (geneGenealogyNameMap.containsKey(egg.getName())) {
	// 	throw (new RuntimeException("ERROR: duplicate gene genealogy name " + egg.getName() + ". Check inputs and try again."));
	//     }
	    
	//     geneGenealogyNameMap.put(egg.getName(), geneGenealogy);
	// }

    protected void readSwitchingFrequency (BufferedReader inbr) throws Exception {
	System.out.println("Across-row switching frequency gamma: ");
	double inGamma = Double.parseDouble(inbr.readLine().trim());
	gamma = new SwitchingFrequency(SwitchingFrequency.GAMMA, inGamma, calculationCache, parentalTreeNameMap.keys().size());
    }

    /**
     * Read transition probability parameters (other than those related to basic coalescent model calculations, i.e., parental tree
     * branch lengths).
     *
     * Current model uses two parameters: a hybridization frequency $v$, and a recombination frequency $u$.
     * See writeup for details.
     *
     * Need to change this into a file.
     *
     * To support ratio-solving quickly, need to look-up from (source parental-tree/gene-genealogy names, sink parental-tree/gene-genealogy names) pair to
     * corresponding ratio parameter quickly. Add a map for this.
     */
    protected void readSwitchingFrequencyRatioTermFile (BufferedReader inbr) throws Exception {
	System.out.println("Hidden state switching frequency ratio term input file: ");
	String switchingFrequencyRatioTermFilename = inbr.readLine().trim();

	hiddenStatePairToSwitchingFrequencyRatioTermMap = new BidirectionalMultimap<Tuple<HiddenState,HiddenState>,SwitchingFrequencyRatioTerm>();
	nameToSwitchingFrequencyRatioTermMap = new BijectiveHashtable<String,SwitchingFrequencyRatioTerm>();
	switchingFrequencyRatioTermNameToOptimizeFlagMap = new Hashtable<String,Tuple3<Double,Double,Boolean>>();

	// for now, not used
	// cap on non-self-transition frequencies is part of runHmm.calculateSwitchingFrequencies(...) 
	//
	// to cap switching frequencies
	// so that self-transition frequency never gets too small
	// otherwise model exhibits degenerate behavior
	int numAlternatives = hiddenStates.size();

	BufferedReader br = new BufferedReader(new FileReader(switchingFrequencyRatioTermFilename));
	String line = "";
	while ((line = br.readLine()) != null) {
	    StringTokenizer st = new StringTokenizer(line);
	    if (st.countTokens() < 6) {
		throw (new IOException("ERROR: incorrect number of fields in input line from file " + switchingFrequencyRatioTermFilename + ": " + line));
	    }
	    String name = st.nextToken();
	    Double minimumWeight = new Double(st.nextToken());
	    double initialWeight = Double.parseDouble(st.nextToken());
	    Double maximumWeight = new Double(st.nextToken());
	    Boolean optimizeFlag = new Boolean(st.nextToken());
	    SwitchingFrequencyRatioTerm sfrt = new SwitchingFrequencyRatioTerm(name, initialWeight, calculationCache, numAlternatives);
	    if (nameToSwitchingFrequencyRatioTermMap.containsKey(name)) {
		throw (new IOException("ERROR: duplicate switching frequency ratio term in file " + switchingFrequencyRatioTermFilename + ": " + name));
	    }
	    nameToSwitchingFrequencyRatioTermMap.put(name, sfrt);
	    switchingFrequencyRatioTermNameToOptimizeFlagMap.put(name, new Tuple3<Double,Double,Boolean>(minimumWeight, maximumWeight, optimizeFlag));

	    while (st.hasMoreTokens()) {
		String pairString = st.nextToken();
		StringTokenizer pairTok = new StringTokenizer(pairString, SWITCHING_FREQUENCY_RATIO_TERM_TRANSITION_NAMES_DELIMITER_1);
		if (pairTok.countTokens() != 2) {
		    throw (new IOException("ERROR: incorrect number of hidden states in transition from file " + switchingFrequencyRatioTermFilename + ": " + pairString));
		}
		StringTokenizer hTok1 = new StringTokenizer(pairTok.nextToken(), SWITCHING_FREQUENCY_RATIO_TERM_TRANSITION_NAMES_DELIMITER_2);
		if (hTok1.countTokens() != 2) {
		    throw (new IOException("ERROR: incorrect number of tree names in hidden state from file " + switchingFrequencyRatioTermFilename + ": " + pairString));
		}
		String px = hTok1.nextToken();
		String gx = hTok1.nextToken();
		StringTokenizer hTok2 = new StringTokenizer(pairTok.nextToken(), SWITCHING_FREQUENCY_RATIO_TERM_TRANSITION_NAMES_DELIMITER_2);
		if (hTok2.countTokens() != 2) {
		    throw (new IOException("ERROR: incorrect number of tree names in hidden state from file " + switchingFrequencyRatioTermFilename + ": " + pairString));
		}
		String py = hTok2.nextToken();
		String gy = hTok2.nextToken();
		
		HiddenState hx = parentalTreeGeneGenealogyNamePairToHiddenStateMap.get(px, gx);
		HiddenState hy = parentalTreeGeneGenealogyNamePairToHiddenStateMap.get(py, gy);

		Tuple<HiddenState,HiddenState> pair = new Tuple<HiddenState,HiddenState>(hx, hy);
		// bleh fix this
		if (hiddenStatePairToSwitchingFrequencyRatioTermMap.containsKey(pair) && (hiddenStatePairToSwitchingFrequencyRatioTermMap.get(pair) != null) && (hiddenStatePairToSwitchingFrequencyRatioTermMap.get(pair).size() > 0)) {
		    throw (new IOException("ERROR: duplicate transition for a switching frequency ratio term in file " + switchingFrequencyRatioTermFilename + ": " + sfrt.getName() + " " + pair.toString()));
		}
		hiddenStatePairToSwitchingFrequencyRatioTermMap.put(pair, sfrt);
	    }
	}
    }

    /**
     * Strict!
     * Only call after both switching frequency ratio terms and trees have been built.
     */
    protected boolean verifySwitchingFrequencyRatioTerms () {
	// With new, simpler model,
	// only require terms for within-row switching frequencies.
	//
	// All possible ordered pairs of hidden states must be specified, 
	// *except* for ordered pairs where both members are the same parental tree 
	// (corresponding to transitions
	// where the parental tree remains unchanged).
	for (int i = 0; i < hiddenStates.size(); i++) {
	    HiddenState hi = hiddenStates.get(i);
	    for (int j = 0; j < hiddenStates.size(); j++) {
		HiddenState hj = hiddenStates.get(j);
		
		// comment this out later if adding 
		// across-row switching frequency parameterization
		if (!checkSameParentalTreeClass(hi, hj)) {
		    continue;
		}
		
		Tuple<HiddenState,HiddenState> pair = new Tuple<HiddenState,HiddenState>(hi, hj);
		if (!hiddenStatePairToSwitchingFrequencyRatioTermMap.containsKey(pair)) {
		    System.err.println ("ERROR: hiddenStatePairToSwitchingFrequencyRatioTermMap is missing an entry for hidden state pair " + pair);
		    return (false);
		}
	    }
	}

	return (true);
    }

//		if (i != j) {

    protected void readSubstitutionModelParameters (BufferedReader br) throws Exception {
	System.out.println("Input non-zero substitution model rates in format <AG> <AC> <AT> <GC> <GT>: ");
	String line = br.readLine();
	StringTokenizer st = new StringTokenizer(line);
	// // Let GTRSubstitutionModel take care of checking number of input rates.
	// if (st.countTokens() != GTRSubstitutionModel.getSubstitutionParameterCount()) {
	//     throw (new IOException("ERROR: invalid number of substitution model rates."));
	// }
	double[] rates = new double[st.countTokens()];
	int i = 0;
	while (st.hasMoreTokens()) {
	    rates[i] = Double.parseDouble(st.nextToken());
	    i++;
	}
	
	for (double rate : rates) {
	    if (rate <= 0.0) {
		throw (new IOException("ERROR: substitution model rates must be strictly positive."));
	    }
	}

	System.out.println("Input non-zero substitution model base frequencies in format <A> <G> <C> <T>: ");
	line = br.readLine();
	st = new StringTokenizer(line);
	if (st.countTokens() != NucleotideAlphabet.getClassInstance().getAlphabet().length()) {
	    throw (new IOException("ERROR: invalid number of substitution model base frequencies."));
	}
	double[] baseFrequencies = new double[st.countTokens()];
	i = 0;
	while (st.hasMoreTokens()) {
	    baseFrequencies[i] = Double.parseDouble(st.nextToken());
	    i++;
	}

	for (double baseFrequency : baseFrequencies) {
	    if (baseFrequency <= 0.0) {
		throw (new IOException("ERROR: substitution model base frequencies must be strictly positive."));
	    }
	}

	gtrSubstitutionModel = new GTRSubstitutionModel();
	gtrSubstitutionModel.setSubstitutionRates(rates, baseFrequencies);
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

		// //check += a[i][j];
		// //System.out.println ("inner loop in calculateAij(): " + a[i][j]);
		// if (checkSameParentalTreeClass(si, sj)) {
		//     a[i][j] *= (1.0 - transitionProbabilityParameters.getHybridizationFrequency());
		// }
		// else {
		//     // need to divide by the number of parental tree classes other than the current parental tree class
		//     int numParentalTreeClasses = parentalTreeClasses.keySet().size();
		//     // shouldn't come through in this case, 
		//     // but paranoid
		//     if (numParentalTreeClasses <= 1) {
		// 	throw (new RuntimeException("ERROR: different parental classes in calculateAij() but number of parental classes is less than two."));
		//     }
		    
		//     a[i][j] *= transitionProbabilityParameters.getHybridizationFrequency() / (numParentalTreeClasses - 1);
		// }

    // kliu - no longer needed
    /**
     * Given the current model parameter values, 
     * calculate the maximum possible value for a frequency parameter.
     * Anything greater will cause a self-transition probability to go negative.
     *
     * Hmm... need to expose this externally.
     */
    // public double calculateMaximumFrequencyParameter (TransitionProbabilityParameters.ParameterChoice parameterChoice) {
    // 	double minmax = -1.0;
    // 	boolean initializedFlag = false; // bleh

    // 	// strict!
    // 	if (hiddenStates.size() <= 1) {
    // 	    throw (new RuntimeException("ERROR: not enough hidden states in runHmm.calculateMaximumFrequencyParameter(...)."));
    // 	}

    // 	// take min max value over all states
    // 	for (int i = 0; i < hiddenStates.size(); i++) {
    // 	    HiddenState si = hiddenStates.get(i);
    // 	    double sumCoalescentContributionRecombination = 0.0;
    // 	    double sumCoalescentContributionHybridization = 0.0;

    // 	    // compute max for a state si
    // 	    for (int j = 0 ; j < hiddenStates.size(); j++) {
    // 		// set self-transition probability at the end
    // 		if (i == j) {
    // 		    continue;
    // 		}
		
    // 		HiddenState sj = hiddenStates.get(j);
    // 		if (checkSameParentalClass(si, sj)) {
    // 		    sumCoalescentContributionRecombination += sj.calculateProbabilityOfRootedGeneGenealogyInParentalTree();
    // 		}
    // 		else {
    // 		    sumCoalescentContributionHybridization += sj.calculateProbabilityOfRootedGeneGenealogyInParentalTree();
    // 		}
    // 	    }

    // 	    double max;
    // 	    switch (parameterChoice) {
    // 	    case RECOMBINATION_FREQUENCY:
    // 		max = (1 - sumCoalescentContributionHybridization * transitionProbabilityParameters.get(TransitionProbabilityParameters.ParameterChoice.HYBRIDIZATION_FREQUENCY)) / sumCoalescentContributionRecombination;
    // 		break;
    // 	    case HYBRIDIZATION_FREQUENCY:
    // 	    default:
    // 		max = (1 - sumCoalescentContributionRecombination * transitionProbabilityParameters.get(TransitionProbabilityParameters.ParameterChoice.RECOMBINATION_FREQUENCY)) / sumCoalescentContributionHybridization;
    // 		break;
    // 	    }

    // 	    // update appropriately
    // 	    if (!initializedFlag || (max < minmax)) {
    // 		minmax = max;
    // 		initializedFlag = true;
    // 	    }
    // 	}

    // 	return (minmax);
    // }

