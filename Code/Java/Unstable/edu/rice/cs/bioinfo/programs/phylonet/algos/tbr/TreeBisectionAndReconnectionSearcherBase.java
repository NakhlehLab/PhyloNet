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

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.UnrootedTreespaceHeuristicSearcherBase;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/23/12
 * Time: 9:58 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TreeBisectionAndReconnectionSearcherBase<T, N, E, S> extends UnrootedTreespaceHeuristicSearcherBase<T, N, E, S>
    {
        protected /*final*/ Func2<T, E, Boolean> isInternalEdge;

        protected /*final*/ Func2<T, N, Boolean> isLeaf;

        public TreeBisectionAndReconnectionSearcherBase(Predicate1<T> isRooted, Func3<T, N, E, Boolean> isDestinationNode, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, final Func2<T, N, Iterable<E>> getIncidentEdges, final Func2<T, E, Tuple<N, N>> getNodesOfEdge)
        {
        	super(isRooted, isDestinationNode, getNodes, getEdges, getIncidentEdges, getNodesOfEdge);

            isLeaf = new Func2<T, N, Boolean>() {
                public Boolean execute(T tree, N node) {
                    return  IterableHelp.countInt(getIncidentEdges.execute(tree, node)) == 1;
                }
            };

            isInternalEdge = new Func2<T, E, Boolean>() {
                public Boolean execute(T tree, E edge) {
                    Tuple<N,N> edgeNodes = getNodesOfEdge.execute(tree, edge);
                    return !isLeaf.execute(tree, edgeNodes.Item1) && !isLeaf.execute(tree, edgeNodes.Item2);
                }
            };

        }

        protected /*override*/ abstract T searchWithoutTreeValidation(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations, S givenTreeScore);

        public /*override*/ void assertValidTree(T tree, Predicate1<T> isRooted, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N, N>> getNodesOfEdge)
        {
            super.assertValidTree(tree, isRooted, getNodes, getEdges, getIncidentEdges, getNodesOfEdge);

            for (N node : getNodes.execute(tree))
        {
            Iterator<E> incidentEdges = getIncidentEdges.execute(tree, node).iterator();
            incidentEdges.next();
            if (incidentEdges.hasNext())
            {
                return;
            }
        }

        throw new IllegalArgumentException("Given tree does not contain an internal edge.");
        }
    }