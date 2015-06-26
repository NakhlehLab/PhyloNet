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

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.StringWriter;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * This class implements the methodes declared in the interface Network.
 *
 * @author Cuong Than
 *
 * @param <T> indicates the type of additional data this node will store.
 *
 * June 12, 06: Created
 */
public class BniNetwork<T> implements Network<T>, Cloneable {
	/**
	 * This constructor creates an empty network.
	 */
	public BniNetwork()
	{
		_root = null;
	}

	/**
	 * This constructor creates a network with <code>root</code> as the network's root.
	 *
	 * @param root: The node that will becomes the root of the network to be created.
	 */
	public BniNetwork(BniNetNode<T> root)
	{
		if (root != null && root.getIndeg() > 0) {
			System.err.println("The root of the network has predecessors.");
			_root = null;
		}
		else {
			_root = root;
		}

	}

	/**
	 * This function creates the root for a network.
	 *
	 * @return <code>true</code> if the network is initially empty. In this case,
	 * a new root for the network will be created. Otherwise, the fuction returns <code>false</code>.
	 */
	public boolean createRoot(String name)
	{
		if (isEmpty()) {
			_root = new BniNetNode<T>();
			_root.setName(name);
			return true;
		}
		else {
			return false;
		}
	}

	/**
	 * @return the (unique) root of the network.
	 */
	public NetNode<T> getRoot()
	{
		return _root;
	}

	/**
	 * Return the node in the network with <code>name</code>. If no such node
	 * exists, return null.
	 */
	public NetNode<T> findNode(String name)
	{
		for (NetNode<T> node : bfs()) {
			if (node.getName().equals(name)) {
				return node;
			}
		}

		return null;
	}

	/**
	 * @return <code>true</code> if the network is empty; <code>false</code> otherwise.
	 */
	public boolean isEmpty()
	{
		return (_root == null);
	}

	/**
	 * @return a BFS iterable list of all nodes in the network.
	 */
	public Iterable<NetNode<T>> bfs()
	{
		return new BfsSearch<T>(_root);
	}

	/**
	 * @return a DFS iterable list of all nodes in the network.
	 */
	public Iterable<NetNode<T>> dfs()
	{
		return new DfsSearch<T>(_root);
	}

	/**
	 * @return an iterable list of leaves of the network.
	 */
	public Iterable<NetNode<T>> getLeaves()
	{
		List<NetNode<T>> leaves = new LinkedList<NetNode<T>>();
		for (NetNode<T> node : bfs()) {
			if (node.isLeaf()) {
				leaves.add(node);
			}
		}

		return leaves;
	}

	/**
	 * @return a set of network nodes.
	 */
	public Iterable<NetNode<T>> getNetworkNodes()
	{
		LinkedList<NetNode<T>> networkNodes = new LinkedList<NetNode<T>>();
		for (NetNode<T> node : dfs()) {
			if (node.isNetworkNode()) {
				networkNodes.add(node);
			}
		}

		return networkNodes;
	}

	/**
	 * @return a set of tree nodes.
	 */
	public Iterable<NetNode<T>> getTreeNodes()
	{
		LinkedList<NetNode<T>> treeNodes = new LinkedList<NetNode<T>>();
		for (NetNode<T> node : bfs()) {
			if (node.isTreeNode()) {
				treeNodes.add(node);
			}
		}

		return treeNodes;
	}

	/**
	 * This function clears all nodes in the network.
	 */
	public void clear()
	{
		_root = null;
	}

	/**
	 * This function creates an identical copy of this current network.
	 */
    public Network<T> clone()
    {
        /*
        try {
            return (Network) super.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException("This should be impossible ...");
        }
        */

        return Networks.readNetwork(this.toString());
    }
    /*
	public Network<T> clone()
	{
		BniNetNode<T> root_copy = new BniNetNode<T>();

		return new BniNetwork<T>(root_copy);
	}
	*/

	/**
	 * This function checks if the network has duplicate names or not.
	 *
	 * @return: <code>true</code> if the network contains duplicate names;
	 * <code>false</code> otherwise.
	 */
	public boolean hasDuplicateNames()
	{
		Map<String, NetNode<T>> map = new Hashtable<String, NetNode<T>>();

		for (NetNode<T> node : bfs()) {
			String name = node.getName();

			if (name != NetNode.NO_NAME) {
				if (map.get(name) == null) {
					map.put(name, node);	// First instance of such name.
				}
				else {
					return true;			// Seconde instance of that name.
				}
			}
		}

		return false;	// All node names are distinct.
	}


    /**
     * @return the number of leaves in the network.
     */
    public int getLeafCount(){
        int count = 0;
        for(NetNode node: this.getLeaves()){
            count ++;
        }
        return count;
    }


    /**
     * @return the number of reticulations in the network.
     */
    public int getReticulationCount(){
        int count = 0;
        for(NetNode node: this.getNetworkNodes()){
            count ++;
        }
        return count;
    }


    /**
     * @return the number of edges in the network.
     */
    public int getEdgeCount(){
        int count = 0;
        for(NetNode node: Networks.postTraversal(this)){
            count += node.getChildCount();
        }
        return count;
    }


    /**
     * Reset the root of the network
     */
    public void resetRoot(NetNode root){
        _root = (BniNetNode)root;
    }


    public String toString(){
        RnNewickPrinter<T> rnNewickPrinter = new RnNewickPrinter<T>();
        StringWriter sw = new StringWriter();
        rnNewickPrinter.print(this, sw);
        return sw.toString();
    }


	// Data members
	private BniNetNode<T> _root;
}
