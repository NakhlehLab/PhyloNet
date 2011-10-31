package edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.*;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

/**
 * This class allows us to visit all nodes under a start node in a network.
 * All the nodes, including the start node, are visited in breadth-first search.
 * 
 * @author Cuong Than
 */
public class BfsSearch<T> implements Iterable<NetNode<T>> {
	// Constructors
	public BfsSearch(NetNode<T> start)
	{
		_start_node = start;
	}
	
	public Iterator<NetNode<T>> iterator()
	{
		return new BfsIterator();
	}

	// Data memebers
	private NetNode<T> _start_node;
	
	// Iterator helper class
	private class BfsIterator implements Iterator<NetNode<T>> {
		// Constructors
		public BfsIterator()
		{
			_queue = new LinkedList<NetNode<T>>();				
			_marked_nodes = new LinkedList<NetNode<T>>();
			_current_node = null;
		
			// Initialize the queue.
			if (_start_node != null) {				
				if (!_queue.offer(_start_node)) {
					System.err.println("Cannot initialize the BFS queue.");
					return;
				}
			}
		}
		
		/**
		 * This function checks if there are any unvisited nodes.
		 */
		public boolean hasNext()
		{
			return !_queue.isEmpty();
		}
		
		/**
		 * This function returns the next element in the network.
		 */
		public NetNode<T> next()
		{
			assert(!_queue.isEmpty());				
			
			// Update the queue and the list of marked nodes.
			_current_node = _queue.poll();				
			for (NetNode<T> node : _current_node.getChildren()) {
				if (node.isTreeNode()) {
					_queue.offer(node);
				}
				else if (!_marked_nodes.contains(node)) {
					_queue.offer(node);
					_marked_nodes.add(node);	// Only need to prevent network nodes from multiple visits.
				}
			}
			
			return _current_node;
		}
		
		/**
		 * This method removes a node from network.
		 */
		public void remove()
		{			
			System.err.println("This method is currently not supported.");
			return;
		}
		
		// Data members
		private Queue<NetNode<T>> _queue;		// Stores nodes to be visited.
		private List<NetNode<T>> _marked_nodes;	// Stores marked nodes.
		private NetNode<T> _current_node;		// The currently visited node.
	}
}