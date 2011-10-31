package edu.rice.cs.bioinfo.programs.phylonet.structs.network.util;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/6/11
 * Time: 6:13 PM
 * To change this template use File | Settings | File Templates.
 */

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;

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
