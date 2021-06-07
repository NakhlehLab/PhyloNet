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

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;

import java.io.Serializable;

/**
 * This interface defines methods that a network must implement.
 *
 * @author Cuong Than and Derek Ruths
 *
 * May 27, 06: Created
 */
public interface Network<T> extends Serializable {
	/**
	 * This function creates the root for a network.
	 *
	 * @return <code>true</code> if the network is initially empty. In this case,
	 * a new root for the network will be created. Otherwise, the fuction returns <code>false</code>.
	 */
	public boolean createRoot(String name);

	/**
	 * @return the (unique) root of the network. It returns <code>null<code> if the network
	 * is empty.
	 */
	public NetNode<T> getRoot();

	/**
	 * @return <code>true</code> if the network is empty; <code>false</code> otherwise.
	 */
	public boolean isEmpty();

	/**
	 * Return the node in the network with <code>name</code>.
	 */
	public NetNode<T> findNode(String name);

	/**
	 * @return an iterable list of all leaves in the network.
	 */
	public Iterable<NetNode<T>> getLeaves();

	/**
	 * @return a BFS iterable list of all nodes in the network.
	 */
	public Iterable<NetNode<T>> bfs();

	/**
	 * @return a DFS iterable list of all nodes in the network.
	 */
	public Iterable<NetNode<T>> dfs();

	/**
	 * @return a set of network nodes, i.e. list of nodes with more than one parent.
	 */
	public Iterable<NetNode<T>> getNetworkNodes();

	/**
	 * @return a set of tree nodes, i.e. list of nodes with exactly one parent.
	 */
	public Iterable<NetNode<T>> getTreeNodes();

	/**
	 * This function clears all nodes in the network.
	 */
	public void clear();

	/**
	 * This function creates an identical copy of this current network.
	 */
	public Network<T> clone();

	/**
	 * This function checks if the network has duplicate names or not.
	 *
	 * @return: <code>true</code> if the network contains duplicate names;
	 * <code>false</code> otherwise.
	 */
	public boolean hasDuplicateNames();


    /**
     * @return the number of leaves in the network.
     */
    public int getLeafCount();


    /**
     * @return the number of reticulations in the network.
     */
    public int getReticulationCount();

    /**
     * @return the number of edges in the network.
     */
    public int getEdgeCount();

    /**
     * Reset the root of the network
     */
    public void resetRoot(NetNode root);

    /**
     * @return the Rich Newick string of the network.
     */
    public String toString();

}

