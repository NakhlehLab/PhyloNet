package edu.rice.bioinfo.programs.phylonet.structs.tree.model;

import java.util.List;

public interface TNode {

	public static final double NO_DISTANCE = Double.NEGATIVE_INFINITY;

	public static final String NO_NAME = "";

	/**
	 * @return <code>true</code> if this node is an ancestor of <code>node</code>.
	 */
	public boolean isAncestor(TNode node);

	/**
	 * @return <code>true</code> if this node is presently a member of a tree.  This method
	 * returns <code>false</code> after a node has been removed from the tree.
	 */
	public boolean isValid();

	/**
	 * @return the unique ID of this node.  This ID cannot be changed during the lifetime of this tree.
	 */
	public int getID();

	/**
	 * @return the name of this node
	 */
	public String getName();

	/**
	 * Return the tree to which this node belongs.  If the node has been removed
	 * from the tree, then this method will return <code>null</code>.
	 */
	public Tree getTree();

	/**
	 * @return the parent of this node.  Only the root may have
	 * a <code>null</code> parent.
	 */
	public TNode getParent();

	/**
	 * @return the distance between this node and its parent.
	 */
	public double getParentDistance();

	/**
	 * set the distance between this node and its parent.
	 */
	public void setParentDistance(double distance);

	/**
	 * @return this node's children.
	 */
	public Iterable<? extends TNode> getChildren();

	/**
	 * @return the number of children this node has.
	 */
	public int getChildCount();

	/**
	 * @return <code>true</code> if this node is a leaf node and has no children.
	 */
	public boolean isLeaf();

	/**
	 * @return <code>true</code> if this node is the root of the tree.
	 */
	public boolean isRoot();

	/**
	 * @return the number of leaves below this node.
	 */
	public int getLeafCount();

	public List<TNode> getSiblings();

	/**
	 *
	 * @return an iterable list of leaves below this node.
	 */
	public Iterable<TNode> getLeaves();

	/**
	 *
	 * @return an iterable list of nodes below this node (including this node) that are visited
	 * in the post order fashion.
	 */
	public Iterable<TNode> postTraverse();
}