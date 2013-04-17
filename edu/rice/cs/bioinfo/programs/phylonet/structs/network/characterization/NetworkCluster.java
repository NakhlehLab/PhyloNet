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

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.BitSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/28/11
 * Time: 3:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkCluster<T>
{
    /**
	 * This function constructs a cluster for a given node in a tree generated from a network.
	 *
	 * @param net: The network this cluster belongs to.
	 * @param nodes: The set of leaves in the cluster.
	 */
	public NetworkCluster(Network<T> net, List<NetNode<T>> nodes)
	{
		_net = net;
		init(nodes);
	}

	/**
	 * This function initializes the cluster by converting nodes into a BitSet.
	 */
	private void init(List<NetNode<T>> nodes)
	{
		// Initialize the bit vector holding nodes in this cluster.
		_bits = new BitSet();
		_bits_size = 0;

		for (NetNode<T> leaf : _net.getLeaves()) {
			if (nodes != null && nodes.contains(leaf)) {
				_bits.set(_bits_size);
			}
			_bits_size++;
		}

		// Store the size of this cluster.
		if (nodes != null) {
			_size = nodes.size();
		}
		else {
			_size = 0;
		}
	}

	/**
	 * This function returns all the leaves in this cluster.
	 */
	public List<NetNode<T>> getLeaves()
	{
		List<NetNode<T>> nodes = new LinkedList<NetNode<T>>();
		int i = 0;

		for (NetNode<T> leaf : _net.getLeaves()) {
			if (_bits.get(i)) {
				nodes.add(leaf);
			}
			i++;
		}

		return nodes;
	}

	/**
	 * This function returns the names of all leaves in this cluster.
	 */
	public List<String> getNames()
	{
		List<String> names = new LinkedList<String>();
		int i = 0;

		for (NetNode<T> leaf : _net.getLeaves()) {
			if (_bits.get(i)) {
				names.add(leaf.getName());
			}
			i++;
		}

		return names;
	}

	/**
	 * This method combines leaves of this cluster and <code>nc</code>. This function assumes
	 * that two clusters have the same list of leaves, both in order and the leaves themselves.
	 *
	 * @param nc: The other NetworkCluster to be combined with this cluster.
	 */
	public void union(NetworkCluster<T> nc)
	{
		_bits.or(nc._bits);	// Update the bits representing nodes in the cluster.
		_size += nc._size;	// Update its size.
	}

	/**
	 * This method checks if two clusters are equal or not. This function assumes
	 * that two clusters have the same list of leaves, both in order and the leaves themselves.
	 *
	 * @param nc: The other cluster to be compared.
	 *
	 * @return: <code>true</code> if the two networks are equal; <code>false</code> otherwise.
	 */
	public boolean equals(Object nc)
	{
		assert(nc instanceof NetworkCluster);

		boolean identical;
		NetworkCluster<T> ref = (NetworkCluster<T>) nc;

		if (_net != ref._net) {	// Two clusters are from two networks.
			List<String> names = getNames();
			List<String> refnames = ref.getNames();

			if (names.size() != refnames.size()) {
				identical = false;
			}
			else {
				identical = true;
				for (String str : names) {
					if (!refnames.contains(str)) {
						identical = false;
						break;
					}
				}
			}
		}
		else {
			identical = _bits.equals(ref._bits);
		}

		return identical;
	}

    public int hashCode(){
        int code = 0;
        for(String name: getNames()){
            code += name.hashCode();
        }
        return code;
    }

	/**
	 * This function tests if this cluster is empty, i.e. there're no leaves in it.
	 */
	public boolean isEmpty()
	{
		return (_size == 0);
	}

	/**
	 * This function tests if a cluster is trivial. A cluster is trivial if it contains
	 * just one leaf or it contains all the leaves in the network.
	 */
	public boolean isTrivial()
	{
		return (_size == 1 || _size == _bits_size);
	}

	/**
	 * This function prints this cluster.
	 */
	public String toString()
	{
		Iterator<NetNode<T>> it = getLeaves().iterator();
		String str = new String();

		while (it.hasNext()) {
			str += it.next().getName();
			if (it.hasNext()) {
				str += " ";
			}
		}

		return str;
	}

	// Data members
	private Network<T> _net;	// The network containing this cluster.
	private BitSet _bits;		// The BitSet corresponding to leaves in a cluster.
	private int _size;			// Holds the number of leaves in this cluster.
	private int _bits_size;		// Holds the number of leaves in the tree.
}
