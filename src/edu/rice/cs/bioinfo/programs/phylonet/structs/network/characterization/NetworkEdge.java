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

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/28/11
 * Time: 3:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkEdge<T>
{
    // Constructors
	public NetworkEdge(NetNode<T> tail, NetNode<T> head)
	{
		_tail = tail;
		_head = head;
	}

	/**
	 * This function returns the tail of this edge.
	 */
	public NetNode<T> getTail()
	{
		return _tail;
	}

	/**
	 * This function returns the head of this edge.
	 */
	public NetNode<T> getHead()
	{
		return _head;
	}

	/**
	 * This function tests if the two edges are equal.
	 *
	 * @param ne: The edge to be tested.
	 */
	public boolean equals(Object ne)
	{
		assert(ne instanceof NetworkEdge);

		NetworkEdge ref = (NetworkEdge) ne;

		return (_tail == ref._tail && _head == ref._head);
	}

	// Data members
	private NetNode<T> _head;	// The head of this edge.
	private NetNode<T> _tail;	// The tail, i.e. parent of _head.
}
