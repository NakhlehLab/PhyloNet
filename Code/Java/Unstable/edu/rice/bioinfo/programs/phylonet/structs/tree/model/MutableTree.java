package edu.rice.bioinfo.programs.phylonet.structs.tree.model;

/**
 * This interface defines the methods that a tree that can be changed.
 * 
 * @author Derek Ruths
 * 
 * TODO: Add tree change listeners
 */
public interface MutableTree extends Tree {
	
	/**
	 * Set the name of the tree.
	 */
	public void setName(String name);
	
	public void constrainByLeaves(Iterable<String> leaf_names);
	
	/**
	 * In an empty tree, this method creates a new node which is the root.
	 * 
	 * @return the new root.
	 * 
	 * @throws RuntimeException if the tree is not empty.
	 */
	public TMutableNode createRoot();

	// overridden Tree methods
	public TMutableNode getRoot();
	
	public TMutableNode getNode(int id);
	
	public TMutableNode getNode(String name);
	
	public Iterable<? extends TMutableNode> getNodes();
	
	public String toStringWD();
}
