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

package edu.rice.cs.bioinfo.programs.phylonet.algos.nni;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.uci.ics.jung.graph.Graph;

import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 4/10/12
 * Time: 7:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class NNISearcherNeighborsExpansiveInPlaceJUNG<N,E,S>
{
     Func1<Graph<N,E>,Iterable<N>> getNodes = new Func1<Graph<N,E>,Iterable<N>>()
        {
            public Iterable<N> execute(Graph<N, E> labelGraph) {
                return labelGraph.getVertices();
            }
        };

        Predicate1<Graph<N,E>> isTreeRooted = new Predicate1<Graph<N, E>>() {
            public boolean execute(Graph<N, E> labelGraph) {
                return false;
            }
        };

        Func1<Graph<N,E>,Collection<E>> getEdges = new Func1<Graph<N, E>, Collection<E>>() {
            public Collection<E> execute(Graph<N, E> labelGraph) {
                return labelGraph.getEdges();
            }
        };

        Func2<Graph<N,E>,N,Collection<E>> getIncidentEdges = new Func2<Graph<N,E>,N, Collection<E>>() {
            public Collection<E> execute(Graph<N,E> tree, N label) {
                return tree.getIncidentEdges(label);
            }
        };

}
