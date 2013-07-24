package phylogeny;

import java.util.ArrayList;

// kliu - need to work in log space everywhere

/**
 * Nodes class for a tree
 */
public class Node {

    protected ArrayList<Node> children;
    protected Node parent;
    // kliu - this will need to be set via E-M later
    protected double tbranch;
    protected String taxa;
    protected String obs;
    protected ArrayList<Double> likelihood = new ArrayList<Double>(4);	// likelihood ArrayList
    protected ArrayList<Double> leftLike;
    protected ArrayList<Double> rightLike;
    protected char[] genes = {'A', 'C', 'G', 'T'};
    protected double lamda = .1;		// this value of the evolutionary constant
	
    /**
     * Constructor for Node. 
     * Should use constructors for specific node types instead. LeafNode or InternalNode
     */
    public Node() {
	this.children = null;
	this.parent = null;
	this.tbranch = -1;
	for (int i = 0; i < 4; i++) {
	    likelihood.add(-1.0);
	}
    }
	
    /**
     * Constructor used for cloning an Internal node
     * @param leftChild - left child Node
     * @param rightChild - right child Node
     * @param parent - The parent node
     * @param tbranch - a double that is the branch length from this node to its parent node
     */
    public Node(Node leftChild, Node rightChild, Node parent, double tbranch) {
	this.setLeftChild(leftChild);
	this.setRightChild(rightChild);
	this.parent = parent;
	this.tbranch = tbranch;
	this.taxa = null;
	this.obs = null;
		
	for (int i = 0; i < 4; i++) {
	    likelihood.add(-1.0);
	}
    }
	
    /**
     * Constructor used for cloning a leaf node
     * @param parent - the parent node
     * @param tbranch - a double that is the branch length from this node to its parent node
     */
    public Node(Node parent, double tbranch) {
	this.children = null;
	this.parent = parent;
	this.tbranch = tbranch;
	this.taxa = null;
	this.obs = null;
	for (int i = 0; i < 4; i++) {
	    likelihood.add(-1.0);
	}
    }
	
    /**
     * Constructor for Leaf nodes
     */
    public Node(String taxa, Node parent, double tbranch) {
	this.children = null;
	this.parent = parent;
	this.tbranch = tbranch;
	this.taxa = taxa;
	this.obs = null;
	for (int i = 0; i < 4; i++) {
	    likelihood.add(-1.0);
	}
		
    }
	
    /**
     * Constructor for Leaf Nodes
     */
    public Node(String taxa) {
	this.children = null;
	this.parent = null;
	this.tbranch = -1;
	this.taxa = taxa;
	this.obs = null;
	for (int i = 0; i < 4; i++) {
	    likelihood.add(-1.0);
	}
    }
	
	
    /**
     * Set the Parent node
     * Parent will be null if node is a root
     */
    public void setParent(Node parent) {
	this.parent = parent;
    }
	
    /**
     * Set Left child
     * @param left A Node
     */
    public void setLeftChild(Node left) {
	if (children == null) children = new ArrayList<Node>(2);
	children.add(0, left);
    }
	
    /**
     * Set Right child
     * @param right A Node
     */
    public void setRightChild(Node right) {
	if (children == null) children = new ArrayList<Node>(2);
	children.add(1, right);
    }
	
    /** 
     * Add Child
     * @param child A Node
     */
    public void addChild(Node child) {
	if (children == null) children = new ArrayList<Node>(2);
	children.add(child);
    }
	
    /**
     * Add 2 children
     * @param left - Left child Node
     * @parm right - right child Node
     */
    public void addChildren(Node left, Node right) {
	this.setLeftChild(left);
	this.setRightChild(right);
    }
	
    /**
     * Set time branch from itself to its parent
     * t-branch for root will be a default of -1
     * @param t A Double
     */
    public void setTbranch(double t) {
	this.tbranch = t;
    }

    // kliu - doesn't seem to be any guards about leaf property?
    /**
     * Sets the observation for leaf node
     * @param obs - the observation String {A,T,G,C}
     */
    public void setObs(String obs) {
	this.obs = obs;
		
	// sets likelihood of leaf based on observation
	int index = -1;
	for (int i = 0; i < genes.length; i++) {
	    if (obs.charAt(0) == genes[i]) {
		index = i;
	    }
	}
		
	for (int i = 0; i < 4; i++) {
	    // kliu - these appear to be probabilities, not log probabilities
	    if (i == index) {
		likelihood.set(i, 1.0);
	    }
	    else likelihood.set(i, 0.0);
	}
    }
	
    /**
     * Sets the taxa or species name for the leaf node
     * @param taxa - the taxa String { a species}
     */
    public void setTaxa(String taxa) {
	this.taxa = taxa;
    }
	
    /** 
     * @return Parent Node
     */
    public Node getParent() {
	return parent;
    }
	
    /**
     * @return The Left child node
     */
    public Node getLeft() {
	if (isLeaf())
	    return null;
	else return children.get(0);
    }
	
    /**
     * @return The Right child node
     */
    public Node getRight() {
	if (isLeaf())
	    return null;
	else return children.get(1);
    }
	
    /**
     * @return String Observation {A,T,G,C}
     */
    public String getObs() {
	return obs;
    }
	
    /**
     * @return String Taxa {species}
     */
    public String getTaxa() {
	return taxa;
    }
	
    /** 
     * @return double branch length from node to its parent
     */
    public double getTbranch() {
	return tbranch;
    }
	
	
	
    /**
     * @return The node's likelihood arraylist
     */
    public ArrayList<Double>  getLikelihood() {
	if (isLeaf()) {
	    //System.out.println(" I am leaf : " + taxa + " and my likelihood array is : " + likelihood);
	    return this.likelihood;
	}
	else {
			
	    leftLike = children.get(0).getLikelihood();
	    rightLike = children.get(1).getLikelihood();
			
	    for (int i = 0; i < 4; i++) {
		double tempLeft = 0.0;
		double tempRight = 0.0;
		for (int j = 0; j < 4; j++) {
		    tempLeft += leftLike.get(j) * getPij(genes[i], genes[j], lamda, children.get(0).getTbranch());
		    tempRight += rightLike.get(j) * getPij(genes[i], genes[j], lamda, children.get(1).getTbranch());
		}
		likelihood.set(i, tempLeft * tempRight);
	    }
						
	     //System.out.println(" I am an internal node " + " and my likelihood array is : " + likelihood);

	    return this.likelihood;
	}
    }
	
	
	
	
    /**
     * @return True if node is the root
     */
    public boolean isRoot() {
	if (parent == null) return true;
	return false;
    }
	
	
	
	
    /**
     * @return True if node is a leafNode, false otherwise
     */
    public boolean isLeaf() {
	if (children == null) return true;
	return false;
		
    }
	


    /**
     * Clears the observation and likelihood if is a leaf
     * Otherwise, will just clear any likelihood calculation
     */
    public void clearNode() {
	this.obs = null;
	// kliu - bug here
	//this.likelihood = null;
	// just zero out observations
	for (int i = 0 ; i < likelihood.size(); i++) {
	    likelihood.set(i, 0.0);
	}

	//Recursive step
	if (!this.isLeaf()) {
	    this.getLeft().clearNode();
	    this.getRight().clearNode();
	}
    }
	
	
	

    /**
     * @return An ArrayList of all the LeafNodes that this node can reach
     */
    public ArrayList<Node> getLeaves() {
	if (isLeaf()) {
	    ArrayList<Node> leaf = new ArrayList<Node>();
	    leaf.add(this);
	    return leaf;
	}
	else {
	    ArrayList<Node> leaves = new ArrayList<Node>();
	    leaves.addAll(this.children.get(0).getLeaves());
	    leaves.addAll(this.children.get(1).getLeaves());
	    return leaves;
	}
    }
	
	
	
	
    /**
     * @return A clone of the current node itself
     */
    public Node cloneNode() {
		
	if (isLeaf()) {
	    return new Node(parent.cloneNode(), tbranch);
	}
	else if (isRoot()) {
	    return new Node(children.get(0).cloneNode(), children.get(1).cloneNode(), null, tbranch);
	}
	else {
	    return new Node(children.get(0).cloneNode(), children.get(1).cloneNode(), parent.cloneNode(), tbranch);
	}
    }
	
	

    // kliu - Jukes-Cantor model calculation
    /**
     * @param i - A Nucleotide Letter of type char {A, T, G, C}
     * @param j	- A Nucleotide Letter of type char {A, T, G, C}
     * @param u - A Double - base rate substitution/evolutionary constat
     * @param t - A Double - the time interval or branch length time
     * @return returns the Pij value, or the transition probability between two nucleotides given the branch length
     */
    private double getPij(char i, char j, double u, double t) {
	double p;
        if (!(Character.toString(i)).equals((Character.toString(j)))) {
            return (0.25 + 0.75 * Math.exp((-4.0 / 3.0) * u * t));
        }
        else
            return (0.25 + 0.75 * Math.exp((-4.0 / 3.0) * u * t));
		

    }
	
	
	
	
    /**
     * String representation of the node for debugging and user-viewing purposes
     */
    @Override
    public String toString() {
	if (isLeaf()) {
	    if (obs != null) 
		return taxa + "->" + obs + " :" + tbranch;
	    else
		return taxa + " :" + tbranch;
	}
	else {
	    String toprint = "(" +  children.get(0) + ", " + children.get(1) + ")";
	    if (parent != null) {
		toprint += " :" + tbranch;
	    }
	    return toprint;
	}
    }
	
	
}
