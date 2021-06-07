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

package edu.rice.cs.bioinfo.programs.phylonet.algos.spr;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.UnrootedTreespaceHeuristicSearcherBase;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/23/12
 * Time: 3:04 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SubtreePruningAndRegraftingSearcherBase<T, N, E, S> extends UnrootedTreespaceHeuristicSearcherBase<T, N, E, S>
{
    public SubtreePruningAndRegraftingSearcherBase(Predicate1<T> isRooted, Func3<T, N, E, Boolean> isDestinationNode, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N, N>> getNodesOfEdge)
    {
    	super(isRooted, isDestinationNode, getNodes, getEdges, getIncidentEdges, getNodesOfEdge);

    }

    protected /*override*/ abstract T searchWithoutTreeValidation(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations, S givenTreeScore);

    public /*override*/ void assertValidTree(T tree, Predicate1<T> isRooted, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N, N>> getNodesOfEdge)
    {
        super.assertValidTree(tree, isRooted, getNodes, getEdges, getIncidentEdges, getNodesOfEdge);

        if (IterableHelp.countInt(getEdges.execute(tree)) < 4)
        {
            throw new IllegalArgumentException("Given tree does not contain at least four edges.");
        }
    }
}