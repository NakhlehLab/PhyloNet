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

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.util;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/6/11
 * Time: 6:13 PM
 * To change this template use File | Settings | File Templates.
 */

/**
 * This method implements utility methods on network nodes.
 */
public class NetNodes {
	/**
	 * This function inserts a binary nodes between nodes <code>tail</code> and <code>head</code>.
	 *
	 * @param tail, head: The tail and head node of an edge.
	 * @param tail_distance: The distance from the tail to the new node.
	 *
	 * @return: The node inserted if successfully; null otherwise.
	 */
	public static <T> NetNode<T> breakEdge(NetNode<T> tail, NetNode<T> head, double tail_distance)
	{
		double distance = head.getParentDistance(tail);

		if (distance != NetNode.NO_DISTANCE && distance < tail_distance) {
			System.err.println("Cannot insert a new node.");
			return null;
		}
		else {
			NetNode<T> node = new BniNetNode<T>();

			tail.removeChild(head);
			tail.adoptChild(node, tail_distance);
			if (distance == NetNode.NO_DISTANCE || tail_distance == NetNode.NO_DISTANCE) {
				node.adoptChild(head, NetNode.NO_DISTANCE);
			}
			else {
				node.adoptChild(head, distance - tail_distance);
			}

			return node;
		}
	}
}
