package phylogeny;

import java.util.ArrayList;

/**
 * EvoTree is a class that stores an entire evolution tree.
 * It can be traversed through by recursion - tree/similar to linked list style.
 */
public class EvoTree {
	
	private Node root;					// Tree Root
	private int treeID;					// unique Identification of tree, -1 means it has not yet been identified
	private String aname;
	private int calcLeaves; 			// will be 0 if leaves have not been found, 1 if they have
	private ArrayList<Node> leaves; // stores all the leafNodes of this tree
	
	/**
	 * Constructor for EvoTree
	 * @param root A Node type -- root of the tree
	 */
	public EvoTree(Node root) {
		this.root = root;
		this.treeID = -1;
		this.leaves = null;
		this.calcLeaves = 0;
		this.aname = "";
	}
	
	/**
	 * Empty constructor for EvoTree
	 */
	public EvoTree() {
		this.root = null;
		this.treeID = -1;
		this.leaves = null;
		this.calcLeaves = 0;
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
		this.leaves = null;
		this.calcLeaves = 0; 
	}
	
	/** 
	 * Sets the root of the tree
	 * @param rootIn A Node type
	 */
	public void setRoot(Node rootIn) {
		this.root = rootIn;
		
		// reset leaves --> new root means need to find leaves again
		this.calcLeaves = 0;
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
	
	/**
	 * Returns the probability of this tree given the observed the sequence of genome
	 * P(D|Tree)
	 * @return the Likelihood of this tree with the given leaf sequence
	 * Note: the leaves of the tree must already be populated with the appropriate sequence
	 */
	public double getLikelihood(double u) {
		double result = 0;
		for (int i = 0; i < 4; i++) {
			result += 0.25 * root.getLikelihood(u).get(i);
		}
		return result;
	}
	
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
		// leaves have been found previously already
		if ((calcLeaves == 1) && (leaves != null)) {
			return this.leaves;
		}
		
		// leaves have not yet been found
		else {
			leaves = root.getLeaves();
			calcLeaves = 1;
			return this.leaves;
		}
	}
	
	/**
	 * Creates a string representable version of this tree for the purpose of debugging and reading
	 */
	@Override
	public String toString() {
		return "Tree: " + aname + "\n" + root.toString() + "\n\n";
	}

}
