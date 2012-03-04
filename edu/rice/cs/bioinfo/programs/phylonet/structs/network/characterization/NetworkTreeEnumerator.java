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

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/31/11
 * Time: 3:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkTreeEnumerator<T> implements Iterable<NetworkTree<T>>
{
    	// Constructors
	public NetworkTreeEnumerator(Network<T> net)
	{
		_net = net;
		_network_nodes = new LinkedList<NetNode<T>>();

		for (NetNode<T> node : _net.bfs()) {
			if (node.isNetworkNode()) {
				_network_nodes.add(node);
			}
		}
	}

	public Iterator<NetworkTree<T>> iterator()
	{
		return new NetworkTreeIterator();
	}

	// Data members
	Network<T> _net;	// The network to generate different trees.
	List<NetNode<T>> _network_nodes;	// The set of network nodes of _net.

	// Iterator helper class
	private class NetworkTreeIterator implements Iterator<NetworkTree<T>> {
		// Constructors
		public NetworkTreeIterator()
		{
			_configs = new Stack<Iterator<NetNode<T>>>();
			_edges = new LinkedList<NetworkEdge<T>>();

			// Initialize the first set of network edges chosen to build a tree.
			for (NetNode<T> node : _network_nodes) {
				Iterator<NetNode<T>> it = node.getParents().iterator();
				NetworkEdge<T> ne = new NetworkEdge<T>(it.next(), node);

				_configs.push(it);
				_edges.add(ne);
			}

			assert((_network_nodes.size() == _edges.size()) && (_network_nodes.size() == _configs.size()));
			_has_next = true;
		}

		public boolean hasNext()
		{
			return _has_next;
		}

		public NetworkTree<T> next()
		{
			assert(_has_next);

			// Build the tree from the current configuration.
			List<NetworkEdge<T>> temp = new LinkedList<NetworkEdge<T>>();
			for (NetworkEdge<T> ne : _edges) {
				temp.add(ne);
			}
			NetworkTree<T> tree = new NetworkTree<T>(_net, temp);

			// Prepare the next configuration to build a new tree.
			boolean stop = false;
			int last = _network_nodes.size() - 1;	// Points to the current network node used to build a new config.

			while (!stop && !_configs.isEmpty()) {
				Iterator<NetNode<T>> it = _configs.pop();
				_edges.remove(last--);

				if (it.hasNext()) {
					// Proceed to the next edge in the list of edges from the last removed node to its parents.
					NetNode<T> child = _network_nodes.get(++last);
					_edges.add(new NetworkEdge<T>(it.next(), child));
					_configs.push(it);

					// Reset all edges for nodes after 'last'.
					for (int i = last + 1; i < _network_nodes.size(); i++) {
						child = _network_nodes.get(i);
						it = child.getParents().iterator();

						_edges.add(new NetworkEdge<T>(it.next(), child));
						_configs.push(it);
					}

					stop = true;	// Already found a new configuration.
				}
			}

			// Determine if all trees are already enumerated.
			if (_configs.isEmpty()) {
				_has_next = false;
			}

			return tree;
		}

		public void remove()
		{
			System.err.println("This method is currently not supported.");
			return;
		}

		// Data members
		private Stack<Iterator<NetNode<T>>> _configs;	// This stack is for enumerating all sets of network edges.
		private List<NetworkEdge<T>> _edges;			// List of network edges chosen to build a tree.
		private boolean _has_next;						// True if there're more sets of network edges.
	}
}
