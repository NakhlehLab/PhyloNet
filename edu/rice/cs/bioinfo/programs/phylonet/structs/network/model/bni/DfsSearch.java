package edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.*;

import java.util.*;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.*;

/**
 * This class allows us to visit all nodes under a start node in a network.
 * All the nodes, including the start node, are visited in depth-first search.
 * 
 * @author Cuong Than
 */
public class DfsSearch<T> implements Iterable<NetNode<T>> {
	// Constructors
	public DfsSearch(NetNode<T> start)
	{
		_start_node = start;
	}
	
	public Iterator<NetNode<T>> iterator()
	{
		return new DfsIterator();
	}
	
	// Data members
	NetNode<T> _start_node;
	
	// Iterator helper class
	private class DfsIterator implements Iterator<NetNode<T>> {
		// Constructors
		public DfsIterator()
		{
			_stack = new Stack<NetNode<T>>();
			_marked_nodes = new LinkedList<NetNode<T>>();
			_current_node = null;
			
			// Initialize the stack.
			if (_start_node != null) {									
				_stack.push(_start_node);
			}
		}
		
		/**
		 * This method checks if there are any undiscovered nodes.
		 */
		public boolean hasNext()
		{
			return !_stack.isEmpty();
		}
		
		/**
		 * This method returns the next element in the network.
		 */
		public NetNode<T> next()
		{
			assert(!_stack.isEmpty());
							
			// Update the stack and the list of marked nodes.
			_current_node = _stack.pop();			
			for (NetNode<T> node : _current_node.getChildren()) {
				if (node.isTreeNode()) {
					_stack.push(node);
				}
				else if (!_marked_nodes.contains(node)) {
					_stack.push(node);
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
		private Stack<NetNode<T>> _stack;		// Stores nodes to be visited.
		private List<NetNode<T>> _marked_nodes;	// Stores marked nodes.
		private NetNode<T> _current_node;		// The currently visited node.
	}	
}