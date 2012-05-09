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
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.HashSet;
import java.util.LinkedList;

class TreeValidator
{
    static <T,N,E> void assertValidTree(T tree, Func1<T, Iterable<N>> getNodes, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N,N>> getNodesOfEdge)
    {
        HashSet<N> nodes = new HashSet<N>();
        N someNode = null;

        for (N node : getNodes.execute(tree)) // assert no node duplicates
        {
            nodes.add(node);
            someNode = node;
        }

        if (nodes.size() < 2)
        {
            return;
        }


        HashSet<N> seenNodes = dfsExplore(tree, someNode, getIncidentEdges, getNodesOfEdge);

        if (seenNodes.size() != nodes.size())
        {
            throw new IllegalArgumentException("Given tree is not connected.");
        }


    }

    static <T,N,E> HashSet<N> dfsExplore(T tree, N node, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N,N>> getNodesOfEdge)
    {
        LinkedList<N> nodeNeighbors = new LinkedList<N>();

        for(E edge : getIncidentEdges.execute(tree, node))
        {
            Tuple<N,N> nodesOfEdge = getNodesOfEdge.execute(tree, edge);
            if(nodesOfEdge.Item1.equals(node))
            {
                nodeNeighbors.addLast(nodesOfEdge.Item2);
            }
            else
            {
                nodeNeighbors.addLast(nodesOfEdge.Item1);
            }
        }

        HashSet<N> seenNodes = new HashSet<N>();
        seenNodes.add(node);

        for(N nodeNeighbor : nodeNeighbors)
        {
            TreeValidator.dfsExplore(tree, node, nodeNeighbor, seenNodes, getIncidentEdges, getNodesOfEdge);
        }

        return seenNodes;
    }

    private static <T,N,E> void dfsExplore(T tree, N parent, N child, HashSet<N> seenNodes, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N,N>> getNodesOfEdge)
    {
        if (seenNodes.contains(child))
        {
            throw new IllegalArgumentException("Given tree contains a cycle.");
        }
        seenNodes.add(child);

        LinkedList<N> adjToChildExceptParent = new LinkedList<N>();

        for(E edge : getIncidentEdges.execute(tree, child))
        {
            Tuple<N,N> nodesOfEdge = getNodesOfEdge.execute(tree, edge);
            N otherNode;
            if(nodesOfEdge.Item1.equals(child))
            {
                otherNode =nodesOfEdge.Item2;
            }
            else
            {
                otherNode =nodesOfEdge.Item1;
            }

            if(!otherNode.equals(parent))
            {
                adjToChildExceptParent.addLast(otherNode);
            }
        }


        for(N nonParentChild : adjToChildExceptParent)
        {
            dfsExplore(tree, child, nonParentChild, seenNodes, getIncidentEdges, getNodesOfEdge);
        }
    }
}
