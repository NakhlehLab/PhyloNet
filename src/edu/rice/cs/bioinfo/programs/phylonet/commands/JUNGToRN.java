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

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.uci.ics.jung.graph.Graph;

import java.util.HashSet;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/23/12
 * Time: 5:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class JUNGToRN
{
    public static <N> String toRichNewick(Graph<N,?> tree, boolean isRooted)
    {
        StringBuffer accum = new StringBuffer(";");

        if(tree.getVertexCount() == 0)
        {
            throw new IllegalArgumentException("RichNewick cannot encode an empty tree.");
        }
        else if(tree.getVertexCount() == 1)
        {
            accum.insert(0, tree.getVertices().iterator().next().toString());
        }
        else
        {
            if(isRooted)
            {
               for(N vertex : tree.getVertices())
            {
                if(tree.inDegree(vertex) == 0)
                {
                    HashSet<N> encodedNodes = new HashSet<N>();
                    toRichNewickHelp(vertex, tree, accum, encodedNodes);
                    break;
                }
            }
            }
            else
            {

            for(N vertex : tree.getVertices())
            {
                if(tree.getNeighborCount(vertex) == 3)
                {
                    HashSet<N> encodedNodes = new HashSet<N>();
                    toRichNewickHelp(vertex, tree, accum, encodedNodes);
                    break;
                }
            }
            }
        }


        return (isRooted ? "" : "[&U]") + accum.toString();

    }

    private static <N> void toRichNewickHelp(N vertex, Graph<N, ?> unrootedTree, StringBuffer accum, HashSet<N> encodedNodes) {

        accum.insert(0, vertex.toString());
        encodedNodes.add(vertex);

        HashSet<N> unencodedNeighbors = new HashSet<N>();

        for(N neighbor : unrootedTree.getNeighbors(vertex))
        {
            if(!encodedNodes.contains(neighbor))
            {
                unencodedNeighbors.add(neighbor);
            }
        }

        if(unencodedNeighbors.size() == 0)
        {
            return;
        }
        else
        {
            Iterator<N> neighbors = unencodedNeighbors.iterator();

            accum.insert(0, ")");

            while(neighbors.hasNext())
            {
                N next = neighbors.next();
                toRichNewickHelp(next, unrootedTree, accum, encodedNodes);

                if(neighbors.hasNext())
                {
                    accum.insert(0, ",");
                }
            }
            accum.insert(0, "(");
        }

    }
}
