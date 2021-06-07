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

package edu.rice.cs.bioinfo.programs.phylonet.algos.tbr;

import edu.rice.cs.bioinfo.library.programming.Func1;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/23/12
 * Time: 9:56 PM
 * To change this template use File | Settings | File Templates.
 */
public interface TreeBisectionAndReconstruction<T>
{
    public T performBisection(T tree, Func1<T,Integer> getTreeScore, boolean maximize, int maxIterations);

    public void assertValidTBRTree(T tree);
}
