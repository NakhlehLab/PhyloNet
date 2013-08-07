package phylogeny;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
// kliu - pull in additional library support

// kliu - use this as general tree class.
// Each hidden state corresponds to a (parental tree, gene genealogy) pair.

/**
 * EvoTree is a class that stores an entire evolution tree.
 * It can be traversed through by recursion - tree/similar to linked list style.
 */
public class EvoTree {
    protected static final boolean DEFAULT_DISPLAY_BRANCH_LENGTHS_FLAG = true;
    protected static final boolean DEFAULT_DISPLAY_INTERNAL_NODE_NAMES_FLAG = true;

    protected Node root;					// Tree Root
    protected int treeID;					// unique Identification of tree, -1 means it has not yet been identified
    protected String aname;
    //protected int calcLeaves; 			// will be 0 if leaves have not been found, 1 if they have
    //protected ArrayList<Node> leaves; // stores all the leafNodes of this tree

    /**
     * Constructor for EvoTree
     * @param root A Node type -- root of the tree
     */
    public EvoTree (Node root) {
    this.root = root;
    this.treeID = -1;
    //this.leaves = null;
    //this.calcLeaves = 0;
    this.aname = "";
    }

    /**
     * Empty constructor for EvoTree
     */
    public EvoTree () {
    this.root = null;
    this.treeID = -1;
    //this.leaves = null;
    //this.calcLeaves = 0;
    this.aname = "";
    }

    /**
     * Constructor for EvoTree
     * @param root A Node type -- root of the tree
     * @param treeIDin An Integer that unqiuely defines the tree
     */
    public EvoTree(Node root, int treeIDin) {
    this.root = root;
    this.treeID = treeIDin;
    //this.leaves = null;
    //this.calcLeaves = 0;
    }

    /**
     * Sets the root of the tree
     * @param rootIn A Node type
     */
    public void setRoot(Node rootIn) {
    this.root = rootIn;

    // reset leaves --> new root means need to find leaves again
    //this.calcLeaves = 0;
    }

    /**
     * @return The root of type Node of the tree
     */
    public Node getRoot(){
    return root;
    }

    /**
     * Sets the unique ID for the tree
     * @param idIn An integer that uniquely identifies the tree
     */
    public void setID(int idIn) {
    this.treeID = idIn;
    }

    /**
     * Set name of tree
     */
    public void setName(String name) {
    this.aname = name;
    }

    public String getName() {
    return aname;
    }

    /**
     *
     * @return An integer that uniquely defines this tree
     */
    public int getID() {
    return this.treeID;
    }

    // kliu - entry point for emission probability calculation
    // disable this method - don't want to cache observations on the tree
    // do emission probability calculation on-the-fly
    /**
     * Returns the probability of this tree given the observed the sequence of genome
     * P(D|Tree)
     * @return the Likelihood of this tree with the given leaf sequence
     * Note: the leaves of the tree must already be populated with the appropriate sequence
     */
    // public double getLikelihood() {
    // 	double result = 0;
    // 	for (int i = 0; i < 4; i++) {
    // 	    // kliu - ArrayList recalculated four times over
    // 	    result += 0.25 * root.getLikelihood().get(i);
    // 	}
    // 	return result;
    // }

    /**
     * kliu - On-the-fly version of getLikelihood().
     * Calculates emission probability
     * P[O_t | g(s_i), b_{g(s_i)}, \theta ]
     *
     * See writeup for details.
     */
    public double getLikelihood (ObservationMap column) {
    mapObsToLeaves(column);

    double result = 0;
    ArrayList<Double> rootLikelihoods = root.getLikelihood();
    for (int i = 0; i < 4; i++) {
        // kliu - ArrayList recalculated four times over
        result += 0.25 * rootLikelihoods.get(i);
    }

    // kliu - no need to do this - just overwrite each time
    //clearTree();

    return result;
    }

    // public double getLikelihood (String obs, Map<String, Integer> seqType) {
    // 	mapObsToLeaves(obs, seqType);

    // 	double result = 0;
    // 	for (int i = 0; i < 4; i++) {
    // 	    // kliu - ArrayList recalculated four times over
    // 	    result += 0.25 * root.getLikelihood().get(i);
    // 	}

    // 	clearTree();

    // 	return result;
    // }

    // kliu - gene genealogies associated with second parental
    // tree seem to not have likelihood vector initialized
    // correctly -> NullPointerException
    /**
     * Helper function that Maps observations to leaves.
     */
    protected void mapObsToLeaves (ObservationMap column) {
    ArrayList<Node> leaves = this.getLeaves();
    for (int i = 0; i < leaves.size(); i++) {
        Node leaf = leaves.get(i);
        leaf.setObs(column.get(leaf.getTaxa()));
    }

    // Testing
    //System.out.println("the Tree after mapping obs: " + this);
    }

    // protected void mapObsToLeaves (String obs, Map<String, Integer> seqType) {
    // 	ArrayList<Node> leaves = this.getLeaves();
    // 	for (int i = 0; i < leaves.size(); i++) {
    // 	    Node leaf = leaves.get(i);
    // 	    int index = seqType.get(leaf.getTaxa());
    // 	    String aObs = obs.substring(index,index+1);
    // 	    leaf.setObs(aObs);
    // 	}

    // 	// Testing
    // 	//System.out.println("the Tree after mapping obs: " + this);
    // }

    /**
     * Resets the Tree's leaves (but doesn't reset the leaf types) and the
     * tree nodes' likelihood values
     */
    public void clearTree() {
    root.clearNode();
    }

    /**
     * @return A clone of this current tree
     */
    public EvoTree cloneTree() {
    return new EvoTree(root.cloneNode());
    }

    //	/**
    //	 * @param species A string of length 1 that represents the taxon being added to the tree's leaves
    //	 * @return An arraylist of EvoTrees that are created from adding the new taxon to the old tree
    //	 */
    //	public ArrayList<EvoTree> addNewLeaf(String species) {
    //		// add the new species to this tree
    //		ArrayList<Node> newRoots = root.addNewLeaf(species);
    //
    //		// create new trees from the new roots that resulted from the new addition
    //		ArrayList<EvoTree> newTrees = new ArrayList<EvoTree>();
    //		for (int i = 1; i < newRoots.size(); i++) {
    //			newTrees.add(new EvoTree(newRoots.get(i), i-1));
    //		}
    //		return newTrees;
    //	}

    /**
     *
     * @return An ArrayList of Leaf Nodes of the tree
     */
    public ArrayList<Node> getLeaves() {
    return (root.getLeaves());
    }

    public ArrayList<Node> getNodes() {
    return (root.getNodes());
    }

    /**
     * Get names of taxa.
     */
    public String[] getTaxa () {
    ArrayList<Node> leaves = getLeaves();
    String[] taxa = new String[leaves.size()];
    int i = 0;
    for (Node leaf : leaves) {
        taxa[i] = leaf.getTaxa(); // kliu - change Node.getTaxa() name later
        i++;
    }

    return (taxa);
    }


    // // leaves have been found previously already
    // if ((calcLeaves == 1) && (leaves != null)) {
    //     return this.leaves;
    // }

    // // leaves have not yet been found
    // else {
    //     leaves = root.getLeaves();
    //     calcLeaves = 1;
    //     return this.leaves;
    // }


    /**
     * Creates a string representable version of this tree for the purpose of debugging and reading
     */
    @Override
    public String toString() {
    return "Tree: " + aname + "\n" + root.toString() + "\n\n";
    }

    public String toNewickString () {
    return (toNewickString(DEFAULT_DISPLAY_BRANCH_LENGTHS_FLAG, DEFAULT_DISPLAY_INTERNAL_NODE_NAMES_FLAG));
    }

    /**
     * Stub. Use a flag to indicate whether branch lengths and
     * internal node names should be output or not..
     */
    public String toNewickString (boolean displayBranchLengthsFlag, boolean displayInternalNodeNamesFlag) {
        return root.toNewickString(displayBranchLengthsFlag, displayInternalNodeNamesFlag);
    }

    protected static void test (String filename) {
    try {
        BufferedReader ptreesbr = new BufferedReader(new FileReader(filename));
        TreeParser ptp = new TreeParser(ptreesbr);
        ArrayList<EvoTree> trees = ptp.nexusFileTreeNames(filename);
        ptreesbr.close();

        int i = 0;
        for (EvoTree tree : trees) {
        System.out.println ("Tree " + i + ": |" + tree.toNewickString(true, true) + "|");
        }
    }
    catch (IOException ioe) {
        System.err.println (ioe);
        System.exit(1);
    }
    }

    protected static void test2 (String filename) {
    try {
        BufferedReader ptreesbr = new BufferedReader(new FileReader(filename));
        TreeParser ptp = new TreeParser(ptreesbr);
        ArrayList<EvoTree> trees = ptp.nexusFileTreeNames(filename);
        ptreesbr.close();

        EvoTree tree = trees.get(0);

        if (trees.size() < 1) {
        System.err.println("ERROR: must be at least one tree in input file " + filename + ". Aborting.");
        System.exit(1);
        }
        else if (trees.size() > 1) {
        System.err.println ("WARNING: more than one tree in input file " + filename + ". Only using first tree.");
        }

        for (Node n : tree.getNodes()) {
        System.out.println (n.getTaxa() + " " + n.getTbranch() + " " );
        }
        System.out.println();

        HashMap<String,String> hm = new HashMap<String,String>();
        hm.put("human", "A");
        hm.put("chimp", "A");
        hm.put("gorilla", "G");
        ObservationMap column = new ObservationMap(hm);
        System.out.println ("Column likelihood: |" + tree.getLikelihood(column) + "|");
    }
    catch (IOException ioe) {
        System.err.println (ioe);
        System.exit(1);
    }
    }

    public static void main (String[] args) {
    if (args.length != 1) {
        System.err.println ("Usage: java phylogeny.EvoTree <input tree file in Nexus format>");
        System.exit(1);
    }
    test(args[0]);
    }
}
