package edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model;

/**
 * This class defines the methods and attributes of a node in a {@link MutableTree}.  A TMutableNode belongs
 * to the tree for which it was created.  Any attempt to merge the node with another tree
 * results in a RuntimeException being thrown.  Once a node has been removed from the tree,
 * it becomes <I>invalid</I> for the remainder of its life.
 * 
 * @author Derek Ruths
 */
public interface TMutableNode extends TNode {

	/**
	 * Set the name of this node.
	 * 
	 * @param name is the new name of this node.
	 */
	public void setName(String name);
	
	/**
	 * Set the distance between this node and its parent.
	 */
	public void setParentDistance(double distance);
	
	/**
	 * Create a new child beneath this node.
	 * 
	 * @return the new child.
	 */
	public TMutableNode createChild();
	
	/**
	 * Create a new child with a specific name beneath this node.
	 * 
	 * @param name is the name of the node
	 * 
	 * @return the child created.
	 */
	public TMutableNode createChild(String name);
	
	/**
	 * Create a new child and any subchildren as a child of this node based on
	 * a pre-existing node either in this tree or in some other tree.  A copy of
	 * the node and all nodes beneath it are placed as a child of this node.  The
	 * only property that is guaranteed to be copied is the node name and node
	 * structure. 
	 * 
	 * @param clade is the pre-existing sub-tree.
	 * 
	 * @return the new child.
	 */
	public TMutableNode createChild(TNode clade);
	
	/**
	 * Move an existing node beneath this node.
	 * 
	 * @param nchild is the child for this node to adopt.
	 */
	public void adoptChild(TMutableNode nchild);
	
	/**
	 * This method removes this node and all its children from the tree.
	 */
	public void removeNode();
	
	/**
	 * Remove the specified child from this node.
	 * 
	 * @param child is the node to remove
	 * @param adopt_all if this term is <code>true</code> then this node will adopt all
	 * the children of the child it removes.
	 */
	public void removeChild(TMutableNode child, boolean adopt_all);
	
	/**
	 * Reorganize the tree to which this node belongs so that this node is the root of the tree.
	 */
	public void makeRoot();
	
	// overriden methods from TNode
	public Iterable<? extends TMutableNode> getChildren();
	
	public TMutableNode getParent();
	
	public MutableTree getTree();
}
