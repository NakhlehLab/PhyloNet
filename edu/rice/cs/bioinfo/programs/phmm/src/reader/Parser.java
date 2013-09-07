package reader;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import util.Constants;
import phylogeny.EvoTree;
import phylogeny.Node;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
//import be.ac.ulg.montefiore.run.jahmm.MyHMM;
// kliu - pull in additional library support

public class Parser {
    public static final String TAXON_NAME_LINE_MARKER = ">";

    protected ArrayList<String> alphabet;								/* ArrayList of legal characters/symbols */
    protected ArrayList<ObservationMap> sequence;					/* a REUSABLE arraylist that holds the converted and final usable sequence of integers for ONE file */
    protected ArrayList<HiddenState> trees_states;						/* the list of states in the hmm */
    protected HashMap<String, Integer> myHashmap;						/* mapping from the alphabet to Integers */
    protected HashMap<String, Integer> seqTypes;						/* mapping from types of sequences to integers */
    protected Vector<String> taxa; //  reverse the above map
    protected HashMap<String, String> alleleSpeciesMap;					/* Allele to Species Mapping */
    protected int seqNum;												/* Number of sequences */

    /**
     * Constructor For the Parser
     * Needs the first basic Information file in order to start reading genome sequences
     *
     * @limit Number-of-sequences 10
     * @Assumptions Number-of-symbols less than or equal to 4
     * @param basicFileName - String representation of the file name of the basic file to be parsed
     *
     * <b><u>File format:</u></b> <br />
     * numberOfSequences (--should be same as number of species)<br />
     * Legal Alphabet symbols each separated by a space<br />
     * Species or Taxa each separated by a space (note: should appear in same order as they are in sequence files)<br /><br />
     *
     * <b><u> Example Basic File </u></b> <br />
     * 3
     * A C T G
     * Human Chimp Gorilla
     * @throws IOException
     */

    public Parser(String basicFileName) throws Exception {
	alphabet = new ArrayList<String>();
	myHashmap = new HashMap<String, Integer>();
	seqTypes = new HashMap<String, Integer>();
	taxa = new Vector<String>();
	trees_states = null;
	seqNum = 0;

	String[] line;
	String read;

	// Reads in the necessary initialization information:
	BufferedReader br = new BufferedReader(new FileReader(basicFileName));

	// Read in the first line: number of sequences
	if ((read = br.readLine()) != null) {
	    // Save data as variables
	    seqNum = Integer.parseInt(read);
	} else {
	    br.close();
	    throw new ParserFileException("Error while reading number of Sequences in file.");
	}

	// Read in alphabet symbols
	if ((read = br.readLine()) != null) {
	    line = read.split(" ");
	    for (int i = 0; i < line.length; i++) {
		alphabet.add(line[i]);
	    }
	} else {
	    br.close();
	    throw new ParserFileException("Error while reading alphabet symbols in file.");
	}


	System.out.println("The legal alphabet symbols: " + alphabet);



	//build Hashmap for integer translation from the alphabet
	myHashmap = IntTranslator.buildMap(alphabet);

	// Read in Sequence Types (taxa or species)
	if ((read = br.readLine()) != null) {
	    line = read.split(" ");
	    System.out.println("The number of species: " + seqNum);
	    System.out.print("The species: ");
	    for (int i = 0; i < line.length; i++) {
		System.out.print(line[i] + " ");
		seqTypes.put(line[i], i);
		taxa.add(line[i]);
	    }
	} else {
	    br.close();
	    throw new ParserFileException("Error while reading sequence type or species or taxa in file.");
	}

	// strict!
	if ((taxa.size() != seqNum) || (seqTypes.size() != seqNum)) {
	    br.close();
	    throw new ParserFileException("ERROR: actual and listed number of taxa in " + basicFileName + " differ.");
	}


	br.close();

    }

    protected void addSequence (String name, StringBuffer sequenceBuffer, Hashtable<String,StringBuffer> map, String filename) {
	// strict!
	if (map.contains(name)) {
	    throw (new RuntimeException("ERROR: sequence file " + filename + " contains duplicate taxon " + name + "."));
	}

	map.put(name, sequenceBuffer);
    }

    protected Hashtable<String,StringBuffer> parseFASTAAlignment (String filename) {
	Hashtable<String,StringBuffer> map = new Hashtable<String,StringBuffer>();
	try {
	    BufferedReader br = new BufferedReader(new FileReader(filename));
	    String line = "";
	    String name = "";
	    StringBuffer sequenceBuffer = null;
	    
	    while ((line = br.readLine()) != null) {
		line = line.trim();
		if (line.startsWith(TAXON_NAME_LINE_MARKER)) {
		    // collect accumulated sequence info and store it now
		    if ((!name.equals("")) && (sequenceBuffer != null)) {
			addSequence(name, sequenceBuffer, map, filename);
		    }
		    
		    // update name
		    // sequence name is the rest of the line, omitting whitespace
		    name = line.substring(1, line.length()).trim();
		    
		    // kliu testing
		    //System.out.println ("name: |" + name + "|");

		    // need a new sequence buffer too
		    sequenceBuffer = new StringBuffer();
		}
		else {
		    // any other line is a sequence line - append it to current sequence buffer
		    sequenceBuffer.append(line);
		}
	    }

	    // do it once more for final taxon, since it isn't followed by another TAXON_NAME_LINE_MARKER
	    if ((!name.equals("")) && (sequenceBuffer != null)) {
		addSequence(name, sequenceBuffer, map, filename);
	    }

	    br.close();
	}
	catch (IOException ioe) {
	    // strict!
	    System.err.println (ioe);
	    ioe.printStackTrace();
	    System.exit(1);
	}

	return (map);
    }

    protected void convertFASTAMapToListOfObservationMap (Hashtable<String,StringBuffer> map, boolean keepUninformativeSitesFlag) {
	sequence = new ArrayList<ObservationMap>();
	// strict!
	int alignmentLength = -1; // undefined
	for (String taxon : map.keySet()) {
	    if ((alignmentLength >= 0) && (map.get(taxon).length() != alignmentLength)) {
		throw (new RuntimeException("ERROR: input sequences are not all of the same length. " + alignmentLength + " " + map.get(taxon).length()));
	    }
	    alignmentLength = map.get(taxon).length();
	}

	// testing
	if (Constants.WARNLEVEL > 1) { System.out.println ("Input alignment length: " + alignmentLength); }

	for (int i = 0; i < alignmentLength; i++) {
	    Hashtable<String,Character> column = new Hashtable<String,Character>();
	    for (String taxon : map.keySet()) {
		column.put(taxon, new Character(map.get(taxon).charAt(i)));
	    }
	    ObservationMap omap = new ObservationMap(column);
	    if (!keepUninformativeSitesFlag && !checkUninformativeSite(omap)) {
		sequence.add(omap);
	    }
	}
    }

    protected boolean checkUninformativeSite (ObservationMap omap) {
	Hashtable<Character,Counter> counts = new Hashtable<Character,Counter>();
	for (String key : omap.keySet()) {
	    Character c = omap.get(key);
	    if (!counts.containsKey(c)) {
		counts.put(c, new Counter());
	    }
	    counts.get(c).increment();
	}

	// eh, informative in the sense of parsimony
	// must be at least two letters each with at least two exhibiting taxa
	//
	// AAAA is not parsimony-informative,
	// but MLE on this column prefers "point" star tree with degenerate zero-length branches
	// to all other trees. In this sense, AAAA is likelihood informative.
	// Meh.
	int numDistinguishingStates = 0;
	for (Character c : counts.keySet()) {
	    if (counts.get(c).get() >= 2) {
		numDistinguishingStates++;
	    }
	}
	
	return (numDistinguishingStates >= 2);
    }

    public void parseMe (String filename, boolean keepUninformativeSitesFlag) {
	Hashtable<String,StringBuffer> map = parseFASTAAlignment(filename);
	convertFASTAMapToListOfObservationMap(map, keepUninformativeSitesFlag);
    }

    // kliu - switch to FASTA format
    /**
     * This function parses a sequence file to extract the
     * observation sequence. In addition, it inputs the emission
     * probability for the each seen observation into the HMM.
     *
     * @param filename genome sequence file to be parsed
     * <p> <b> Input File Format for Each Sequence : </b></p>
     * <i>Num_Of_Sequences Length_of_a_Sequence Alphabet<br>
     * Type_of_Sequence Sequence<br>
     * Type_of_Sequence Sequence<br>
     * ...<br>
     * ... </i><br><br><br>
     *
     * <b> Example </b><br>
     * <i>4 3 actg <br>
     * human atc<br>
     * mouse cgg<br>
     * chimp ggg<br>
     * gorilla ttt </i><br><br>
     *
     */
    // public void parseMe (String filename, Hmm<ObservationMap> myhmm) throws Exception {
    //     List<RandomAccessFile> rafList = new ArrayList<RandomAccessFile>();			/* a list of buffered readers for each sequence */
    //     String read;
    //     int ch;
    //     int seqLen;																/* Length of each sequence */
    //     sequence = new ArrayList<ObservationMap>();									/* The sequence of observation read from current file */


    //     BufferedReader br = new BufferedReader(new FileReader(filename));

    //     // Reads in the length of the sequences
    //     if ((read = br.readLine()) != null) {
    //         // Save data as variables
    //         seqLen = Integer.parseInt(read);
    //     }
    //     else {
    //         br.close();
    //         throw new ParserFileException("File is empty.");
    //     }

    //     br.close();

    //     // Given the number of sequences, initialize the corresponding number of
    //     // buffered readers and get them all to the beginning of their assigned sequence
    //     // in order to read all the sequences at the same time.
    //     long fseekpos = 2;
    //     for (int i = 0; i < seqNum; i++) {

    //         // create new buffers and add to buffer list
    //         rafList.add(new RandomAccessFile(filename, "r"));

    //         // set the current buffer
    //         RandomAccessFile raf = rafList.get(i);

    //         // move start of buffer to the appropriate line
    //         raf.seek(fseekpos);

    //         // read in the type of sequence (eg. human, mouse,...etc)
    //         // and map it to the appropriate line number
    //         //String type = "";		// Uncomment to Test
    //         while ((ch = raf.read() ) != -1) {
    // 		if (ch == ' ') break;
    // 		//else type += (char)ch;		// Uncomment to Test
    // 		fseekpos += 1;
    //         }
    //         if (ch == -1) {
    // 		multiFileCloser(rafList, i); // Closes current file and any file opened
    // 		throw new ParserFileException("File Input is incorrectly formatted!");
    //         }

    //         //now RAF is at the correct position
    //         //set fseekpos for next RAF
    //         fseekpos += (long)seqLen + 1; // add the length of sequence and the \n character


    //         // Testing purposes //
    //         //System.out.println("the type is: " + type);


    //     }

    //     System.out.println("Reading Observation Progress :");

    //     // kliu - looks like order of sequences in sequence file must match
    //     // order in basic info file
    //     // get students to change this to FASTA or something
    //     // remove basic info file

    //     // Read the observations column by column
    //     for (int j = 0; j < seqLen; j++) {
    //         //String obs = "";
    //         int letter;
    //         Hashtable<String,Character> column = new Hashtable<String,Character>();
    //         for (int i = 0; i < rafList.size(); i++) {
    //             if ((letter = rafList.get(i).read()) != -1) {
    //                 //obs += (char)letter;
    //                 column.put(taxa.get(i), new Character((char) letter));
    //             }
    // 		else {
    // 		    multiFileCloser(rafList, rafList.size()-1);
    // 		    throw new ParserFileException("Error while reading observation sequence!");
    // 		}
    //         }


    //         //convert to corresponding integer sequence
    //         //int obsInt = IntTranslator.letterToInt(obs, myHashmap);

    //         // kliu - bleh - just keep around all columns


    //         // add Integer obs to sequence
    //         //sequence.add(obsInt);
    //         sequence.add(new ObservationMap(column));

    //         // kliu - no longer necessary - OpdfMap objects created as part of HMM creation in MyHMM.buildMyHMM
    //         //input emissions probabilities for observations into every state/tree
    //         //inputEmissions(obs, obsInt, trees_states, seqTypes, myhmm);

    //         //keep a progress report
    //         if ((j%10000) == 0) {
    // 		System.out.print("\r"+ (((double)j/(double)seqLen) * 100.0) + "%");
    //         }

    //         //			if (j % 1000000 == 0) {
    //         //				MemoryReport.report();
    //         //			}


    //     }


    //     System.out.print("\r 100.0% Done!                     ");
    //     System.out.println("\n");

    //     multiFileCloser(rafList, rafList.size()-1);


    //     /* Testing purposes */
    //     //////////////////////
    //     //		for (int i = 0; i < sequence.size(); i++) {
    //     //			System.out.println("\n " + sequence.get(i));
    //     //			System.out.println("\n" + IntTranslator.intToLetter(seqNum,sequence.get(i), myHashmap));
    //     //		}
    //     /////////////////////

    // }


    /**
     * Function to parse Allele to Species Mapping
     * @param filename - filename path of file to be parse
     */
    public void parseAlleleSpecies (String filename) throws Exception {
        alleleSpeciesMap = new HashMap<String,String>();
        int numLines;
        String[] line;
        BufferedReader br = new BufferedReader(new FileReader(filename));

        //Read in the number of lines
        numLines = Integer.parseInt(br.readLine());

        // Reads in mapping
        for (int i = 0; i < numLines; i++) {
            line = br.readLine().split(" ");
            alleleSpeciesMap.put(line[0], line[1]);
        }

        br.close();
    }

    /**
     * Function to set the list of trees/states of the hmm model
     * @param trees - an arraylist of type EvoTree
     */
    public void setTrees(ArrayList<HiddenState> trees) {
        this.trees_states = trees;
    }

    /**
     * @return the Allele to Species Mapping
     */
    public HashMap<String,String> getAlleleSpeciesMap() {
        return alleleSpeciesMap;
    }

    /**
     * @return Number of Taxa or types of sequences
     */
    public int getNumSeq() {
        return seqNum;
    }

    public List<String> getGeneGenealogyTaxa () {
	return (taxa);
    }

    /**
     * @return List of legal symbols/alphabet
     */
    public ArrayList<String> getAlphabet(){
        return alphabet;
    }

    /**
     * @return Alphabet-Integer mapping
     */
    public HashMap<String,Integer> getAImap() {
        return myHashmap;
    }

    /**
     * @return Types of Genomes/Sequences-Integer mapping
     */
    public HashMap<String, Integer> getTypeMap() {
        return seqTypes;
    }


    /**
     * @return the final converted Observation sequence for use in HMM
     */
    public ArrayList<ObservationMap> getObs() {
        return sequence;
    }


    /**
     *
     * Closes bufferedReaders in the given input list starting from 0th index up to and including
     * the given index.
     * @param brList A list type that contains bufferedreaders to be closed
     * @param index The index of the last element to be closed.
     *
     * @return None
     */
    protected void multiFileCloser(List<RandomAccessFile> rafList, int index) {
        int myindex = index;
        if (index >= rafList.size()) myindex = rafList.size()-1;
        for (int i = 0; i <= myindex; i ++) {
            try {
                rafList.get(i).close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Helper function that Maps observations to leaves in all the possible trees
     */
    // public static void mapObsToLeaves(EvoTree aTree, String obs, HashMap<String, Integer> seqType) {
    //     ArrayList<Node> leaves = aTree.getLeaves();
    //     for (int i = 0; i < leaves.size(); i++) {
    //         Node leaf = leaves.get(i);
    //         int index = seqType.get(leaf.getTaxa());
    //         String aObs = obs.substring(index,index+1);
    //         leaf.setObs(aObs);

    //     }

    //     // Testing
    //     //System.out.println("the Tree after mapping obs: " + aTree);

    // }

    // kliu - need to do this on-the-fly as parameter inference proceeds
    // instead of caching up front
    /**
     * Helper function that maps an observation to every tree
     * and calculates the likelihood/ emission probabilities for that observation in each state
     * and inputs that probability into the Hmm
     */
    // public static void inputEmissions(String obs, int obsInt, ArrayList<HiddenState> inputTreesStates, HashMap<String, Integer> seqType, Hmm<ObservationInteger> myhmm) {
    // 	OpdfInteger tempOpdf = (OpdfInteger) myhmm.getOpdf(0);		// only need to check the first state -->
    // 	// because likelihood calculations are done for every state in one go

    // 	// only calculate the likelihood emission probabilities for an observation
    // 	// if it hasn't already been calculated
    // 	if (tempOpdf.probability(new ObservationInteger(obsInt)) > 1.0) {
    // 	    for (int i = 0; i < inputTreesStates.size(); i++) {
    // 		EvoTree geneTree = inputTreesStates.get(i).getGeneGenealogy();
    // 		//mapObsToLeaves(geneTree, obs, seqType);
    // 		// kliu - emission probability calculation is in here
    // 		// input Likelihood into HMM here
    // 		//double myLikelihood = geneTree.getLikelihood();
    // 		MyHMM.setEmission(myhmm, i, obsInt, geneTree.getLikelihood(obs, seqType));
    // 		// kliu - just clears observation and likelihood cache in EvoTree container of Node objects
    // 		//geneTree.clearTree();
    // 	    }
    // 	}

    // }


    // bleh
    public class Counter {
	protected int i;
	public Counter () {
	    i = 0;
	}

	public void increment () {
	    i++;
	}

	public int get () {
	    return (i);
	}
    }
}
