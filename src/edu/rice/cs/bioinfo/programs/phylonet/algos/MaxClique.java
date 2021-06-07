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

package edu.rice.cs.bioinfo.programs.phylonet.algos;

import java.io.PrintStream;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/3/11
 * Time: 10:54 AM
 * To change this template use File | Settings | File Templates.
 */
public class MaxClique
{
    /**
	 * Used with <code>calculateGroups</code> to indicate that this finder should
	 * find maximal cliques.
	 */
	public static final boolean CLIQUES = false;

	/**
	 * Used with <code>calculateGroups</code> to indicated that this finder should
	 * find maximal independent sets.
	 */
	public static final boolean INDEPENDENT_SETS = true;

	private int matrixSize;

	private boolean currentMode;

	private boolean[][] adjacencyMatrix;

	private int[][] adjacencyLists;

	private int[] numAdjacencies;

	private double[][] data;

	private double currentCutoff;

	private boolean[][] buckets;

	private int[] IS;

	/**
	 * Instantiates a new maximum-clique/independent-set finder using the given
	 * matrix of data.
	 *
	 * @param data	the n-by-n set of distances on which this finder will function
	 */
	public MaxClique(double[][] data) {
		this.data = data;

		matrixSize = data.length;
		if (data[0].length != matrixSize) {
			throw new IllegalArgumentException("non-square data matrix");
		}

		adjacencyMatrix = new boolean[matrixSize][matrixSize];
		numAdjacencies = new int[matrixSize];
		adjacencyLists = new int[matrixSize][matrixSize];
		IS = new int[matrixSize];
		buckets = new boolean[matrixSize][matrixSize];
	}

	/**
	 * Returns the set of maximal cliques or maximal independent sets, determined
	 * by the supplied cutoff value.
	 *
	 * @param mode		<code>MaxClique.CLIQUES</code> if this finder should find
	 * 					cliques, or <code>MaxClique.INDEPENDENT_SETS</code> if
	 * 					this finder should find independent sets
	 * @param cutoff	the value above which nodes are considered adjacent
	 * @return			a list of the indices of the nodes comprising the maximal
	 * 					independent sets or cliques
	 */
	public List<int[]> calculateGroups(boolean mode, double cutoff) {
		currentMode = mode;
		currentCutoff = cutoff;

		calculateAdjacencyMatrix();
		calculateAdjacencyListsAndStuff();

		List<int[]> results = new LinkedList<int[]>();
		backtrack(0, results);
		return results;
	}

	private void calculateAdjacencyMatrix() {
		for(int i=0;i<matrixSize;i++) {
			for(int j=i+1;j<matrixSize;j++) {
				// Items are related if the pairwise value is *greater* than the cutoff.
				adjacencyMatrix[i][j] = (data[i][j] > currentCutoff);
			}
		}
	}

	private void calculateAdjacencyListsAndStuff() {
		for(int i=0;i<matrixSize;i++) {
			for(int j=i+1;j<matrixSize;j++) { // it's i+1, not i, because we ignore self-loops in the graph
				if (adjacencyMatrix[i][j] == currentMode) {
					adjacencyLists[i][numAdjacencies[i]] = j;
					buckets[i][numAdjacencies[i]] = false;
					adjacencyLists[j][numAdjacencies[j]] = i;
					buckets[j][numAdjacencies[j]] = false;
					numAdjacencies[i]++;
					numAdjacencies[j]++;
				}
			}
		}
	}

	private void backtrack(int i, List<int[]> results) {

		if (i >= (matrixSize-1)) {

			// Output new MIS designated by IS
			addIS(results);

		} else {

			int x = i+1;
			int c = 0;

			int y, adjIndexY;
			int z, adjIndexZ;

			// for y \in Adj(x) such that y \leq i
			for (adjIndexY=0; adjIndexY<numAdjacencies[x] && ((y=adjacencyLists[x][adjIndexY]) <= i); adjIndexY++) {
				if (IS[y] == 0) c++;
			}

			if (c == 0) {

				// for y \in Adj(x) such that y \leq i
				for (adjIndexY=0; adjIndexY<numAdjacencies[x] && ((y=adjacencyLists[x][adjIndexY]) <= i); adjIndexY++) {
					IS[y]++;
				}
				backtrack(x, results);
				// for y \in Adj(x) such that y \leq i
				for (adjIndexY=0; adjIndexY<numAdjacencies[x] && ((y=adjacencyLists[x][adjIndexY]) <= i); adjIndexY++) {
					IS[y]--;
				}

			} else {

				IS[x] = c;
				backtrack(x, results);
				IS[x] = 0;
				boolean f = true;
				for (adjIndexY=0; adjIndexY<numAdjacencies[x] && ((y=adjacencyLists[x][adjIndexY]) <= i); adjIndexY++) {
					if (IS[y] == 0) {
						// Put y in Bucket(x):
						buckets[x][adjIndexY] = true;

						// for z \in Adj(y) such that z \leq i
						for (adjIndexZ=0; adjIndexZ<numAdjacencies[y] && ((z=adjacencyLists[y][adjIndexZ]) <= i); adjIndexZ++) {
							IS[z]--;
							if(IS[z] == 0) {
								f = false;
							}
						}
					}
					IS[y]++;
				}
				if (f) {
					backtrack(x, results);
				}
				for (adjIndexY=0; adjIndexY<numAdjacencies[x] && ((y=adjacencyLists[x][adjIndexY]) <= i); adjIndexY++) {
					IS[y]--;
				}
				// for y \in Bucket(x) do
				for (adjIndexY=0; adjIndexY<numAdjacencies[x]; adjIndexY++) {
					if (buckets[x][adjIndexY]) {
						y = adjacencyLists[x][adjIndexY];
						for (adjIndexZ=0; adjIndexZ<numAdjacencies[y] && ((z=adjacencyLists[y][adjIndexZ]) <= i); adjIndexZ++) {
							IS[z]++;
						}

						// delete y from Bucket(x)
						buckets[x][adjIndexY] = false;
					}
				}

			}

		}

	}

	private void addIS(List<int[]> results) {
		List<Integer> thisIS = new LinkedList<Integer>();
		for(int i=0; i < matrixSize; i++) {
			// IS[v] = 0 iff v is in MIS (see paper)
			if (IS[i] == 0) {
				thisIS.add(new Integer(i));
			}
		}
		int[] actualThisIS = new int[thisIS.size()];
		for (int i = 0; i < thisIS.size(); i++) {
			actualThisIS[i] = thisIS.get(i).intValue();
		}
		results.add(actualThisIS);
	}

	/**
	 * Prints the current adjacency matrix to the supplied output stream. For testing
	 * purposes.
	 *
	 * @param ps	the <code>PrintStream</code> to which the matrix will be printed
	 */
	public void printMatrix(PrintStream ps) {
		for(int i=0;i<matrixSize;i++) {
			for(int j=0;j<i;j++) {
				ps.print("  ");
			}
			for(int j=i;j<matrixSize;j++) {
				ps.print(adjacencyMatrix[i][j] + " ");
			}
			ps.println();
		}
	}

	/**
	 * Prints the current adjacency lists to the supplied output stream. For testing
	 * purposes.
	 *
	 * @param ps	the <code>PrintStream</code> to which the lists will be printed
	 */
	public void printAdjacencyLists(PrintStream ps) {
		for(int i=0;i<matrixSize;i++) {
			ps.print(i + " ");
			for(int j=0;j<numAdjacencies[i];j++) {
				ps.print(adjacencyLists[i][j] + " ");
			}
			ps.println();
		}
	}
}
