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
package edu.rice.cs.bioinfo.programs.phylonet.algos.bipartitematching;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/3/11
 * Time: 6:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class HungarianBipartiteMatcher {

	// constants
	public static final double NO_EDGE = Double.NEGATIVE_INFINITY;

	private static final int END_OF_PATH = -2;
	private static final int NO_MATCH = -1;

	// constants used in augmenting path calculation
	private static final int NO_PREDECESSOR = -1;
	private static final int SOURCE = -2;

	// fields
	private int _l_size = -1;
	private int _r_size = -1;
	private double[][] _l_weights = null;
	private double[][] _r_weights = null;

	int[][] _matching = null;
	double _matching_weight = 0;

	// data structures for performing augmented path searches
	double[] _l_shortest_path;
	double[] _r_shortest_path;
	int[] _l_predecessor;
	int[] _r_predecessor;

	int _num_edges = 0;

	// constructors
	/**
	 * This constructor allows the bipartite graph to be manually specified.
	 * This specifies the number of nodes in the bipartite graph.  The method {@link addEdge} is used to
	 * specify the connectivity of the graph.
	 *
	 * @param l_size
	 * @param r_size
	 */
	public HungarianBipartiteMatcher(int l_size, int r_size) {
		setDimensions(l_size, r_size);
	}

	// methods
	public int getLeftCount() {
		return _l_size;
	}

	public int getRightCount() {
		return _r_size;
	}

	public void setDimensions(int l_size, int r_size) {
		_l_size = l_size;
		_r_size = r_size;

		resetWeightMatrix();
	}

	/**
	 * Add edge (i,j) to the graph that the matching will be computed over where
	 * <code>i</code> is the index of the left node, <code>j</code> is the
	 * index of the right node, and <code>weight</code> is the weight associated
	 * with the edge.  If <code>weight</code> is {@link NO_EDGE}, then the edge
	 * is removed from this graph if it currently exists.
	 */
	public void addEdge(int i, int j, double weight) {
		if(weight < 0) {
			throw new RuntimeException("Weight must be >= 0");
		}

		if(i < 0 || i >= _l_size) {
			throw new RuntimeException("Invalid i value");
		}

		if(j < 0 || j >= _r_size) {
			throw new RuntimeException("Invalid j value");
		}

		if(weight != NO_EDGE) {
			if(_l_weights[i][j] == NO_EDGE) {
				_num_edges++;
			}

			_l_weights[i][j] = -weight;
		} else {
			if(_l_weights[i][j] != NO_EDGE) {
				_num_edges--;
			}

			_l_weights[i][j] = NO_EDGE;
		}
	}

	/**
	 * @return the weight of the edge between nodes <code>i</code> and <code>j</code>.  If there
	 * is no edge, then the constant {@link NO_EDGE} is returned.
	 */
	public double getEdgeWeight(int i, int j) {

		double weight = (_l_weights[i][j] > _r_weights[j][i])?_l_weights[i][j]:_r_weights[j][i];

		return weight;
	}

	/**
	 * Setup the weight data structure
	 */
	private void resetWeightMatrix() {
		if(_l_size <= 0) {
			throw new RuntimeException("Invalid L size: " + _l_size);
		}

		if(_r_size <= 0) {
			throw new RuntimeException("Invalid R size: " + _r_size);
		}

		_l_weights = new double[_l_size][_r_size];
		_r_weights = new double[_r_size][_l_size];

		// setup the graph to be empty
		for(int i = 0; i < _l_size; i++) {
			Arrays.fill(_l_weights[i], NO_EDGE);
		}

		for(int i = 0; i < _r_size; i++) {
			Arrays.fill(_r_weights[i], NO_EDGE);
		}

		// setup the augmenting path arrays
		_l_shortest_path = new double[_l_size];
		_r_shortest_path = new double[_r_size];
		_l_predecessor = new int[_l_size];
		_r_predecessor = new int[_r_size];
	}

	/**
	 * @return the weight of the matching found.  This method can only be called after a call to {@link findMatching}.
	 */
	public double getMatchingWeight() {
		if(_matching == null) {
			throw new RuntimeException("getMatchingWeight(): Matching has not been calculated yet");
		}

		return _matching_weight;
	}

	/**
	 * @return the actual matching found by the call to {@link findMatching}.  The return value is an Nx2 matrix in which
	 * each row corresponds to an edge that is part of the matching.  <code>[i][0]</code> is the index of the left endpoint
	 * of the ith edge, <code>[i][1]</code> is the right endpoint of the ith edge.
	 */
	public int[][] getMatching() {
		if(_matching == null) {
			throw new RuntimeException("getMatching(): Matching has not been calculated yet");
		}

		int[][] matching = new int[_matching.length][2];

		for(int i = 0; i < matching.length; i++) {
			matching[i][0] = _matching[i][0];
			matching[i][1] = _matching[i][1];
		}

		return matching;
	}

	/**
	 * Compute the maximum weighted matching over the bipartite graph specified using
	 * {@link addEdge}.
	 */
	public void findMatching() {

		if(_num_edges == 0) {
			_matching = new int[0][0];
			_matching_weight = Double.NEGATIVE_INFINITY;
			return;
		}

		int matching_size = 0;

		int[] l_matches = new int[_l_size];
		int[] r_matches = new int[_r_size];

		// start with an empty match
		Arrays.fill(l_matches,-1);
		Arrays.fill(r_matches,-1);

		// loop until no path was augmented
		int[] path = new int[_l_size + _r_size];
		boolean augmentation_found = findAugmentingPath(path, l_matches, r_matches);

		while(augmentation_found) {

			/*
			// what's the path?
			System.out.print("PATH = ");
			for(int k = 0; k < path.length; k++) {
				System.out.print(path[k] + " ");
			}
			System.out.println();
			*/

			// invert the matching of the path
			int i = 0;
			boolean matched = false;
			while(i < (path.length - 1) && path[i+1] != END_OF_PATH) {

				//System.out.print("(" + path[i] + "," + path[i+1] + ") ");

				if(!matched) {
					matching_size++;

					// if the left node path[i] was involved in a match,
					// undo the prior matching
					// This should never happen...
					if(l_matches[path[i]] != NO_MATCH) {
						// sanity check
						//assert false; FIXME See why this assertion is failing, fix it, and reinstate the assertion!
						/*
						_l_weights[path[i]][l_matches[path[i]]] = -_r_weights[l_matches[path[i]]][path[i]];
						_r_weights[l_matches[path[i]]][path[i]] = NO_EDGE;
						 */
					}

					l_matches[path[i]] = path[i+1];
					r_matches[path[i+1]] = path[i];
					_r_weights[path[i+1]][path[i]] = -_l_weights[path[i]][path[i+1]];
					_l_weights[path[i]][path[i+1]] = NO_EDGE;
				} else {
					matching_size--;

					// This is just a sanity check - this should never happen
					if(l_matches[path[i+1]] != path[i]) {
						throw new RuntimeException("BIG PROBLEM");
					}

					//l_matches[path[i+1]] = NO_MATCH;
					//r_matches[path[i]] = NO_MATCH;
					_l_weights[path[i+1]][path[i]] = -_r_weights[path[i]][path[i+1]];
					_r_weights[path[i]][path[i+1]] = NO_EDGE;
				}

				i++;
				matched = !matched;
			}

			//System.out.println();

            /*
			// print the graph
			for(int x = 0; x < _l_size; x++) {
				for(int y = 0; y < _r_size; y++) {
					System.out.println("(" + x + "," + y + ") = " + _l_weights[x][y] + ", " + _r_weights[y][x]);
				}
			}
			*/

			// find the next augmenting path
			augmentation_found = findAugmentingPath(path, l_matches, r_matches);
		}

		// build the matching
		_matching_weight = 0;
		_matching = new int[matching_size][2];
		int idx = 0;

		for(int i = 0; i < _l_size; i++) {
			if(l_matches[i] != NO_MATCH) {
				_matching[idx][0] = i;
				_matching[idx][1] = l_matches[i];
				_matching_weight += _r_weights[l_matches[i]][i];
				idx++;
			}
		}

		return;
	}

	/**
	 * Use the Bellman-Ford algorithm to find the shortest path (by weight) on the
	 * graph.
	 *
	 * @param path is the array in which the nodes of the path will be stored.  Nodes
	 * with an even index are left nodes, nodes with odd indicies are right nodes.
	 * The path either ends at the end of the path or at the first <code>-1</code>
	 * stored in the array.
	 *
	 * @return <code>true</code> if an augmenting path could be found from a left node
	 * to a right node.
	 */
	private boolean findAugmentingPath(int[] path, int[] l_matches, int[] r_matches) {

		// initialize variables used in the calculation
		Arrays.fill(_l_shortest_path,Double.POSITIVE_INFINITY);
		Arrays.fill(_r_shortest_path,Double.POSITIVE_INFINITY);
		Arrays.fill(_l_predecessor, NO_PREDECESSOR);
		Arrays.fill(_r_predecessor, NO_PREDECESSOR);

		for(int i = 0; i < (_l_size + _r_size - 1); i++) {

			// relax all the source's edges (which are to unmatched left nodes)
			for(int j = 0; j < _l_size; j++) {
				if(l_matches[j] == NO_MATCH) { // && _l_shortest_path[j] > 0) {
					_l_shortest_path[j] = 0;
					_l_predecessor[j] = SOURCE;
				}
			}

			// relax all the left and right node's edges
			//System.out.println("Weight size = " + _weights.length + ", " + _weights[0].length);

			for(int j = 0; j < _l_size; j++) {
				for(int k = 0; k < _r_size; k++) {

					// relax the edge (j,k)
					if(_l_weights[j][k] != NO_EDGE) {
						assert _r_weights[k][j] == NO_EDGE;

						if(_r_shortest_path[k] > _l_shortest_path[j] + _l_weights[j][k]) {
							_r_shortest_path[k] = _l_shortest_path[j] + _l_weights[j][k];
							_r_predecessor[k] = j;
						}
					} else if(_r_weights[k][j] != NO_EDGE) { // relax the edge (k,j)
						assert _l_weights[j][k] == NO_EDGE;

						if(_l_shortest_path[j] > _r_shortest_path[k] + _r_weights[k][j]) {
							_l_shortest_path[j] = _r_shortest_path[k] + _r_weights[k][j];
							_l_predecessor[j] = k;
						}
					}
				}
			}

			/*
			// print the graph
			for(int x = 0; x < _l_size; x++) {
				System.out.println("X " + x + " = " + _l_shortest_path[x]);
			}
			for(int y = 0; y < _r_size; y++) {
				System.out.println("Y " + y + " = " + _r_shortest_path[y]);
			}
			*/
		}

		// find the first unmatched node - this is a candidate for the shortest path endpoint
		int shortest_path_idx;
		for(shortest_path_idx = 0; shortest_path_idx < _r_size && r_matches[shortest_path_idx] != NO_MATCH; shortest_path_idx++);

		// if there isn't an unmatched right node
		if(shortest_path_idx == _r_size) {
			return false;
		}

		for(int i = shortest_path_idx + 1; i < _r_size; i++) {
			if(r_matches[i] == NO_MATCH && _r_shortest_path[i] < _r_shortest_path[shortest_path_idx]) {
				shortest_path_idx = i;
			}
		}

		// if there wasn't a shortest path found
		if(_r_shortest_path[shortest_path_idx] == Double.POSITIVE_INFINITY) {
			return false;
		}

		// reconstruct the shortest path
		int idx = shortest_path_idx;
		int path_idx = 0;
		boolean l_side = false;

		while(idx != SOURCE) {
			path[path_idx] = idx;

			if(l_side) {
				idx = _l_predecessor[idx];
			} else {
				idx = _r_predecessor[idx];
			}

			l_side = !l_side;
			path_idx++;
		}

		if(path_idx < path.length) {
			path[path_idx] = END_OF_PATH;
		}

		// reverse path
		int path_end = path_idx - 1;
		for(int i = 0; i < Math.ceil((double)path_idx / 2.0); i++, path_end--) {
			int tmp = path[path_end];
			path[path_end] = path[i];
			path[i] = tmp;
		}

		return true;
	}

}
