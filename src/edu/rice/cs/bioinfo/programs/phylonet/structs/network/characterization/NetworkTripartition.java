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
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BfsSearch;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/31/11
 * Time: 3:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkTripartition<T>
{
    /**
	 * This function constructs a network tripartition.
	 *
	 * @param node: The node to construct this tripartition.
	 * @param leaves: The set of all leaves in the network.
	 */
	public NetworkTripartition(Network<T> net, NetNode<T> node)
	{
		_net = net;
		_node = node;
		init();
	}

	/**
	 * This function builds the tripartition under _node.
	 */
	private void init()
	{
		_strict = new BitSet();
		_reachable = new BitSet();
		Map<NetNode<T>, Integer> access = new Hashtable<NetNode<T>, Integer>();
		List<NetNode<T>> leaves = new LinkedList<NetNode<T>>();

		for (NetNode<T> leaf : _net.getLeaves()) {
			leaves.add(leaf);
		}

		// Pass 1: Init. the set of reachable leaves and the access map.
		for (NetNode<T> nn : new BfsSearch<T>(_node)) {
			int reach = UNDEFINED;

			if (nn.equals(_node)) {
				reach = UNIQUE;		// There's one way for _node to get to itself.
			}
			else if (nn.isTreeNode()) {
				NetNode<T> parent = nn.getParents().iterator().next();	// nn's sole parent.

				if (access.get(parent).intValue() == UNIQUE) {
					reach = UNIQUE;	// The path from _node to nn, through nn's parent, is also unique.
				}
			}

			access.put(nn, new Integer(reach));

			if (nn.isLeaf()) {
				int i = leaves.indexOf(nn);

				_reachable.set(i);
			}
		}

		// Pass 2: Determine reachability for all "undefined" nodes.
		for (NetNode<T> nn : new BfsSearch<T>(_node)) {
			updateReachability(nn, access);
		}

		// Pass 3: Build the set of leaves that are accessible only thru _node.
		for (int i = 0; i < leaves.size(); i++) {
			if (_reachable.get(i)) {
				NetNode<T> leaf = leaves.get(i);	// A reachable leaf under _node.

				if (access.get(leaf).intValue() == UNIQUE) {
					_strict.set(i);
				}
			}
		}
	}

	/**
	 * This function determines if a node is reachable only through _node or it's reachable throuhg _node,
	 * but not necessarily.
	 *
	 * @param nn: The node to be determined its reachability
	 * @param access: The map storing reachabilities of all nodes under _node.
	 *
	 * @return: UNIQUE if <code>nn</code> is accessible only through _node; NONUNIQUE if it is also accessible
	 * from other nodes besides _node.
	 */
	private int updateReachability(NetNode<T> nn, Map<NetNode<T>, Integer> access)
	{
		int reach = access.get(nn).intValue();

		if (reach == UNDEFINED) {
			reach = UNIQUE;
			for (NetNode<T> parent : nn.getParents()) {
				if (access.get(parent) == null) {
					reach = NONUNIQUE;	// This node can be accessed from a parent that is not under _node.
					break;
				}
				else if (updateReachability(parent, access) == NONUNIQUE) {
					reach = NONUNIQUE;	// One of its parent is not necessarily reachable from _node. Neither is it.
					break;
				}
			}

			// Update the map to avoid re-computation.
			access.remove(nn);
			access.put(nn, new Integer(reach));
		}

		return reach;
	}

	/**
	 * This function tests if two network tripartitions are equal or not. This function assumes that
	 * the set of all network leaves used as a reference in the two tripartitions are exactly the same.
	 *
	 * @return: <code>true</code> if the two tripartitions are equal.
	 */
	public boolean equals(Object ntp)
	{
		assert(ntp instanceof NetworkTripartition);

		NetworkTripartition<T> ref = (NetworkTripartition<T>) ntp;
		boolean identical;

		if (_net == ref._net) {	// Two triparitions are from the same network.
			identical = _reachable.equals(ref._reachable) && _strict.equals(ref._strict);
		}
		else {	// They are from different networks. Compare the leaf names.
			List<String> ref_strict_names = new LinkedList<String>();
			List<String> ref_nonstrict_names = new LinkedList<String>();

			for (NetNode<T> leaf : ref.getStrictSubleaves()) {
				ref_strict_names.add(leaf.getName());
			}
			for (NetNode<T> leaf : ref.getNonstrictSubleaves()) {
				ref_nonstrict_names.add(leaf.getName());
			}

			int strict_count = 0;	// Holds the number of strict subleaves.

			// Test if all strict subleaves are in ref's strict subleaves.
			identical = true;
			for (NetNode<T> leaf : getStrictSubleaves()) {
				strict_count++;
				if (!ref_strict_names.contains(leaf.getName())) {
					identical = false;
					break;
				}
			}

			// Prevent the case that set of strict subleaves is a subset of ref's strict subleaves.
			if (strict_count != ref_strict_names.size()) {
				identical = false;
			}

			// Test if all nonstrict subleaves are in ref's nonstrict subleaves.
			if (identical) {
				int nonstrict_count = 0;	// Holds the number of nonstrict subleaves.

				for (NetNode<T> leaf : getNonstrictSubleaves()) {
					nonstrict_count++;
					if (!ref_nonstrict_names.contains(leaf.getName())) {
						identical = false;
						break;
					}
				}

				if (nonstrict_count != ref_nonstrict_names.size()) {
					identical = false;
				}
			}
		}

		return identical;
	}

	/**
	 * This function returns an iterable list of leaves that are not reachable from the root through _node.
	 */
	public Iterable<NetNode<T>> getUnreachableLeaves()
	{
		List<NetNode<T>> temp = new LinkedList<NetNode<T>>();
		int i = 0;

		for (NetNode<T> leaf : _net.getLeaves()) {
			if (!_reachable.get(i)) {	// Non-reachable
				temp.add(leaf);
			}
			i++;
		}

		return temp;
	}

	/**
	 * This function returns an iterable list of leaves that are only reachable from the root through _node.
	 */
	public Iterable<NetNode<T>> getStrictSubleaves()
	{
		List<NetNode<T>> temp = new LinkedList<NetNode<T>>();
		int i = 0;

		for (NetNode<T> leaf : _net.getLeaves()) {
			if (_strict.get(i)) {	// Must go thru _node.
				temp.add(leaf);
			}
			i++;
		}

		return temp;
	}

	/**
	 * This function returns an iterable list of leaves that are reachable from the root through _node, but
	 * not necessarily through it.
	 */
	public Iterable<NetNode<T>> getNonstrictSubleaves()
	{
		List<NetNode<T>> temp = new LinkedList<NetNode<T>>();
		int i = 0;

		for (NetNode<T> leaf : _net.getLeaves()) {
			if (_reachable.get(i) && !_strict.get(i)) {	// Reachable, but not necessarily thru _node
				temp.add(leaf);
			}
			i++;
		}

		return temp;
	}

	/**
	 * This function returns the head node v in edge e = (u, v) used to compute its tripartitions.
	 */
	public NetNode<T> getTripartitionNode()
	{
		return _node;
	}

	/**
	 * This function prints a network to a string. A tripartition is printed in the form <Ae, Be, Ce>.
	 * For a description of Ae, Be and Ce, see Luay's PhyloNet proposal.
	 */
	public String toString()
	{
		Iterator<NetNode<T>> it;
		String str = new String();

		// The set of strict subleaves.
		it = getStrictSubleaves().iterator();
		while (it.hasNext()) {
			str += it.next().getName();
			if (it.hasNext()) {
				str += " ";
			}
		}
		str += "; ";

		// The set of nonstrict subleaves.
		it = getNonstrictSubleaves().iterator();
		while (it.hasNext()) {
			str += it.next().getName();
			if (it.hasNext()) {
				str += " ";
			}
		}
		str += "; ";

		// The set of unreachable leaves.
		it = getUnreachableLeaves().iterator();
		while (it.hasNext()) {
			str += it.next().getName();
			if (it.hasNext()) {
				str += " ";
			}
		}

		return str;
	}

	// Data members
	private NetNode<T> _node;	// The node corresponding to this tripartition.
	private Network<T> _net;	// The network containing this tripartition.
	private BitSet _strict;		// Indicates leaves that are reachable only thru _node.
	private BitSet _reachable;	// Indicates leaves that are reachable, either only or not necessarily, thru _node.

	private final int UNDEFINED = -1;
	private final int NONUNIQUE = 0;
	private final int UNIQUE = 1;
}
