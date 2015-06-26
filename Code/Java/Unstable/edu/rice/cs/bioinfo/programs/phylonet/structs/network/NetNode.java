/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.structs.network;

import java.io.Serializable;

/**
 * This interface defines basic methods that a network node must implement.
 *
 * @author Cuong Than and Derek Ruths
 *
 * May 27, 06: Created
 */
public interface NetNode<T> extends Serializable {
	/**
	 * @return an iterable list of children of the node invoking this method.
	 */
	public Iterable<NetNode<T>> getChildren();

	/**
	 * @return <code>true</code> if the node is the root of the network.
	 */
	public boolean isRoot();

	/**
	 * @return <code>true</code> if the node is a leaf; <code>false</code> otherwise.
	 */
	public boolean isLeaf();

	/**
	 * @return <code>true</code> if this node is a tree node, i.e. it has only one parent.
	 */
	public boolean isTreeNode();

	/**
	 * @return <code>true</code> if this node is a network node, i.e. it has more than one parent.
	 */
	public boolean isNetworkNode();

	/**
	 * @return an iterable list of parents of the node invoking this method.
	 */
	public Iterable<NetNode<T>> getParents();


    /**
     * @return if parent is a parent node
     */
    public boolean hasParent(NetNode<T> parent);

	/**
	 * @param parent: The parent of this node.
	 *
	 * @return the distance from this node to <code>parent</code>; NaN if <code>parent</code>
	 * is not in the list of parents of this node.
	 */
	public double getParentDistance(NetNode<T> parent);

    /**
     * Returns the probability property of the edge between this node and the given parent.
     *
     * @param parent The parent that defines the in-edge to this node containing the desired probability property.
     * @return  The [0...1] probability property.
     * @throws IllegalArgumentException if the given node is not a parent of the node.
     */
    public double getParentProbability(NetNode<T> parent);

    /**
     * Sets the probability property of the edge between this node and the given parent.
     *
     * @param parent The parent that defines the in-edge to this node containing the desired probability property.
     * @param  probability The new probability value.
     * @throws IllegalArgumentException if the given node is not a parent of the node.
     */
    public void setParentProbability(NetNode<T> parent, double probability);

    /**
     * Returns the support property of the edge between this node and the given parent.
     *
     * @param parent The parent that defines the in-edge to this node containing the desired support property.
     * @return  The [0...1] probability property.
     * @throws IllegalArgumentException if the given node is not a parent of the node.
     */
    public double getParentSupport(NetNode<T> parent);

    /**
     * Sets the support property of the edge between this node and the given parent.
     *
     * @param parent The parent that defines the in-edge to this node containing the desired support property.
     * @param  support The new support value.
     * @throws IllegalArgumentException if the given node is not a parent of the node.
     */
    public void setParentSupport(NetNode<T> parent, double support);

	/**
	 * @return in-degree of a node, which equals to the number of parents of this node.
	 */
	public int getIndeg();

	/**
	 * @return osut-degree of a node, which equals to the number of children of this node.
	 */
	public int getOutdeg();

	/**
	 * @return the name of this node.
	 */
	public String getName();

	/**
	 * @return data stored in this node.
	 */
	public T getData();

	/**
	 * This function connects an existing node (the node that makes a call to this
	 * function) to another node <code>child</code>. The calling code will add <code>child</code>
	 * to its list of children if <code>child</code> has not been already a child of the calling node.
	 *
	 * @param child: The node that the calling node wants to connects to.
	 * @param distance: The distance between the calling code and <code>child</code>.
	 *
	 * @return <code>true</code> if this function succeeded; <code>false</code> otherwise.
	 */
	public boolean adoptChild(NetNode<T> child, double distance);

	/**
	 * This function makes <code>child</code> no longer a child of this node.
	 *
	 * @param child: The node to be removed.
	 *
	 * @return: <code>true</code> if <code>child</code> is indeed a child of this node.
	 */
	public boolean removeChild(NetNode<T> child);


    public void removeItself();


	/**
	 * This function changes the name of an existing node.
	 *
	 * @param name: The new name for the calling node.
	 */
	public void setName(String name);

	/**
	 * This function updates data stored in this node.
	 *
	 * @param data is the new data.
	 */
	public void setData(T data);

	/**
	 * This functions sets the distance from this calling node to <code>parent</code> with
	 * the new value <code>newDistance</code>.
	 *
	 * @param parent: A parent of this node that it wants to modify the distance.
	 * @param distance: New value for the distance from this node to <code>parent</code>.
	 *
	 * @return: <code>true</code> if the operation parent is indeed a parent of this node;
	 * <code>false</code> otherwise.
	 */
	public boolean setParentDistance(NetNode<T> parent, double distance);

	/**
	 * This function returns the number of parent of a node.
	 */
	public int getParentCount();


    /**
     * This function returns the number of children of a node.
     */
    public int getChildCount();


	// Data members
	public static final double NO_DISTANCE = Double.NEGATIVE_INFINITY;	// Constant for no distances.
    public static final double NO_PROBABILITY = Double.NEGATIVE_INFINITY;	// Constant for no probability.
    public static final double NO_SUPPORT = Double.NEGATIVE_INFINITY;// Constant for no support.
	public static final String NO_NAME = "";		// Constant for nodes without a name.
}
