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

package edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;

/**
 * This class models a single HGT event.  An HGT event occurs between two edges in a tree.  This is encoded
 * by providing the child of the edge involved in the HGT event.  Therefore, if the HGT source was the edge
 * between nodes p->c (where p is the parent and c is the child), then the source will be reported as
 * <code>c</code>.
 *
 * In the notation of this object, <code>source</code> is the node which accepts the node and <code>destination</code>
 * is the donator of the node.
 *
 * @author Derek Ruths, Cuong Than
 */
public class HgtEvent implements Serializable {

	// fields
	private TNode _dest;
	private TNode _src;
	private boolean _bad;		// Indicates if an event is bad or not.
	private boolean _violated;	// Indicates if an event is violated or not.

	protected List<TNode> _decomp_nodes;		// List of decomposition nodes this event belong to.
	protected TNode _primary_node;

	// constructors
	public HgtEvent(TNode src, TNode dest) {
		_src = src;
		_dest = dest;

		_decomp_nodes = new LinkedList<TNode>();
		_primary_node = null;

		checkBad();
		checkViolation();
	}

	// methods
	protected void setPrimaryNode(TNode node)
	{
		_primary_node = node;
	}

	protected void addDecompNode(TNode node)
	{
		_decomp_nodes.add(node);
	}

	/**
	 * Return the number of alternatives for this event.
	 */
	protected int getNumberOfAlternatives()
	{
		int count = 1;

		if (_primary_node != null) {
			count *= _primary_node.getParent().getChildCount();
		}

		return count;
	}

	/**
	 * If an event is bad, this function returns <code>true</code>. See comments
	 * for function checkBad for what a bad event is.
	 *
	 * @return <code>true</code> if this event is bad.
	 */
	public boolean isBad()
	{
		return _bad;
	}

	/**
	 * If an event is time-violated, this function returns <code>true</code>. See
	 * comments for function checkViolation for what a violated event is.
	 *
	 * @return <code>true</code> if this event is time-violated.
	 */
	public boolean isViolated()
	{
		return _violated;
	}

	/**
	 * This function checks if an event is bad. An event x -> y is bad if
	 * x is the parent of y. In this case, _bad is set to <code>true</code>.
	 */
	private void checkBad()
	{
		if (_dest.getParent().getName().equals(_src.getName())) {
			_bad = true;
		}
		else {
			_bad = false;
		}
	}

	/**
	 * This function checks if an event is time-violated. An event x -> y is called
	 * time-violated if x is an ancestor, and not parent, of y. In this case,
	 * _violated is set to <code>true</code>.
	 */
	private void checkViolation()
	{
		if (_bad) {
			_violated = false;
		}
		else {
			TNode temp = _dest.getParent();
			_violated = false;

			while (temp != null && !_violated) {
				if (temp.getName().equals(_src.getName())) {
					_violated = true;
				}
				else {
					temp = temp.getParent();
				}
			}
		}
	}

	public TNode getSourceEdge() {
		return _src;
	}

	public TNode getDestEdge() {
		return _dest;
	}

	public int hashCode() {
		return _src.getName().hashCode() + _dest.getName().hashCode();
	}

	public boolean equals(Object obj) {
		if(obj instanceof HgtEvent) {
			HgtEvent event = (HgtEvent) obj;

			return event.getSourceEdge().getName().equals(_src.getName()) && event.getDestEdge().getName().equals(_dest.getName());
		} else {
			return false;
		}
	}

	public String toString()
	{
		String str = _src.getName() + " -> " + _dest.getName() + " ";

		if (_primary_node != null) {
			String name = _primary_node.getName();
			str += "[" + name + "*]";
		}

		for (TNode dnode : _decomp_nodes) {
			String name = dnode.getName();
			str += "[" + name + "]";
		}

		return str;
	}
}

