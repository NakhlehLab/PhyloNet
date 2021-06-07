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

package edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model;

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

	/**
	 * Prune the tree so that it only contains leaves whose names are in the set <code>leaf_names</code>.
	 */
	public void constrainByLeaves(Iterable<String> leaf_names);

	/**
	 * In an empty tree, this method creates a new node which is the root.
	 */
	public TMutableNode createRoot();

	/**
	 * Return the root
	 */
	public TMutableNode getRoot();

	/**
	 * Return the node whose id is given by <code>id</code>
	 */
	public TMutableNode getNode(int id);

	/**
	 * Return the node whose name is given by <code>name</code>
	 */
	public TMutableNode getNode(String name);

	/**
	 * Return all nodes in this tree
	 */
	public Iterable<? extends TMutableNode> getNodes();

	/**
	 * Return the string representation of this tree with its Data
	 */
	public String toStringWD();
}
