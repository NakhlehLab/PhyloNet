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

import edu.rice.cs.bioinfo.library.programming.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/30/12
 * Time: 7:26 PM
 * To change this template use File | Settings | File Templates.
 */
 public abstract class UnrootedTreespaceHeuristicSearcherBase<T, N, E, S> implements UnrootedTreespaceHeuristicSearcher<T,N,E,S>
    {

        protected final Predicate1<T> isRooted;

        protected  final Func3<T, N, E, Boolean> isDestinationNode;

        protected final Func1<T, Iterable<N>> getNodes;

        protected final Func1<T, Iterable<E>> getEdges;

        protected final Func2<T, N, Iterable<E>> getIncidentEdges;

        protected final Func2<T, E, Tuple<N,N>> getNodesOfEdge;

        public UnrootedTreespaceHeuristicSearcherBase(Predicate1<T> isRooted, Func3<T, N, E, Boolean> isDestinationNode, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N,N>> getNodesOfEdge)
        {
            this.isRooted = isRooted;
            this.isDestinationNode = isDestinationNode;
            this.getNodes = getNodes;
            this.getEdges = getEdges;
            this.getIncidentEdges = getIncidentEdges;
            this.getNodesOfEdge = getNodesOfEdge;
        }

        public T search(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations)
        {
            assertValidTree(tree, isRooted, getNodes, getEdges, getIncidentEdges, getNodesOfEdge);
            return searchWithoutTreeValidation(tree, getTreeScore, isBetterScore, maxIterations);
        }


        public T searchWithoutTreeValidation(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations)
        {
            return searchWithoutTreeValidation(tree, getTreeScore, isBetterScore, maxIterations, getTreeScore.execute(tree));
        }

        protected abstract T searchWithoutTreeValidation(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations, S givenTreeScore);

        public /*override*/ void assertValidTree(T tree, Predicate1<T> isRooted, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N, N>> getNodesOfEdge)
        {
            TreeValidator.assertValidTree(tree, getNodes, getIncidentEdges, getNodesOfEdge);
        }
    }

