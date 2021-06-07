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

package edu.rice.cs.bioinfo.programs.phylonet.algos.recomp;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 1/24/12
 * Time: 6:12 PM
 * To change this template use File | Settings | File Templates.
 */
public interface WindowComparator {
    /**
	 * This method compares two adjacent windows.  <code>window1</code> preceeds
	 * <code>window2</code>.  The comparison score does not have to be normalized.
	 *
	 * @param window1_coord is the starting coordinate of window 1
	 * @param window2_coord is the starting coordinate of window 2
	 * @param window1 is the subsequence of each sequence within the first window
	 * @param window2 is the subsequence of each sequence within the second window
	 * @return the comparison score.
	 */
	public double compare(int window1_coord, int window2_coord, char[][] window1, char[][] window2);

	/**
	 * If this method returns <code>true</code>, then the comparison function must be inverted (1-score) after normalization
	 * in order for the high value to be taken as an indicator of recombination.
	 */
	public boolean invertNormalizedScore();
}
