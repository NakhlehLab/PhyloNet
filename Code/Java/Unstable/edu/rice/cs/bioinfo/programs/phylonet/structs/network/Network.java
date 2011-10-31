package edu.rice.cs.bioinfo.programs.phylonet.structs.network;

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
}

