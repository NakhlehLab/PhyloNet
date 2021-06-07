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

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.Tuple;


/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/30/12
 * Time: 7:23 PM
 * To change this template use File | Settings | File Templates.
 */
public interface UnrootedTreespaceHeuristicSearcher<T, N, E, S>
{
        T search(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations);

        T searchWithoutTreeValidation(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations);

        public /*override*/ void assertValidTree(T tree, Predicate1<T> isRooted, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N, N>> getNodesOfEdge);
}
