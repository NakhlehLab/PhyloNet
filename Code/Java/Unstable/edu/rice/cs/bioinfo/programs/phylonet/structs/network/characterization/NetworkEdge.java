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
