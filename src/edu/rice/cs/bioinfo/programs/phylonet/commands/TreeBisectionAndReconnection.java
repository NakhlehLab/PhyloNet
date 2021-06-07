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

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.tbr.TBRSearchNeighborsExpansiveInPlace;
import edu.uci.ics.jung.algorithms.filters.KNeighborhoodFilter;
import edu.uci.ics.jung.graph.Graph;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/23/12
 * Time: 4:54 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("tbreconnection")
public class TreeBisectionAndReconnection extends HeuristicTreeSearchBase
{

     TreeBisectionAndReconnection(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                                  Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader, Random rand) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader, rand);
    }





    @Override
    protected String produceResult() {
         Func1<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,Iterable<NetworkToJUNG.Label>> getNodes = new Func1<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,Iterable<NetworkToJUNG.Label>>()
        {
            public Iterable<NetworkToJUNG.Label> execute(Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]> labelGraph) {
                return labelGraph.getVertices();
            }
        };



        Func3<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>, NetworkToJUNG.Label, NetworkToJUNG.Label[], Boolean> isDestinationNode = new Func3<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>, NetworkToJUNG.Label, NetworkToJUNG.Label[], Boolean>() {
            public Boolean execute(Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> tree, NetworkToJUNG.Label node, NetworkToJUNG.Label[] edge) {
                return edge[1] == node;
            }
        };



        Func1<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,Collection<NetworkToJUNG.Label[]>> getEdges = new Func1<Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]>, Collection<NetworkToJUNG.Label[]>>() {
            public Collection<NetworkToJUNG.Label[]> execute(Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]> labelGraph) {
                return labelGraph.getEdges();
            }
        };

        Func2<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,NetworkToJUNG.Label,Collection<NetworkToJUNG.Label[]>> getIncidentEdges = new Func2<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,NetworkToJUNG.Label, Collection<NetworkToJUNG.Label[]>>() {
            public Collection<NetworkToJUNG.Label[]> execute(Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> tree, NetworkToJUNG.Label label) {
                return tree.getIncidentEdges(label);
            }
        };

        Func2<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,NetworkToJUNG.Label[],Tuple<NetworkToJUNG.Label,NetworkToJUNG.Label>> getNodesOfEdge = new Func2<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,NetworkToJUNG.Label[],Tuple<NetworkToJUNG.Label,NetworkToJUNG.Label>>()
        {
            public Tuple<NetworkToJUNG.Label, NetworkToJUNG.Label> execute(Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> tree, NetworkToJUNG.Label[] edge)
            {
                return new Tuple<NetworkToJUNG.Label, NetworkToJUNG.Label>(edge[0], edge[1]);
            }
        };

        Proc2<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,NetworkToJUNG.Label[]> removeEdge = new Proc2<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,NetworkToJUNG.Label[]>()
        {

            public void execute(Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> tree, NetworkToJUNG.Label[] edge) {

                tree.removeEdge(edge);

            }
        };

        Func3<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,NetworkToJUNG.Label,NetworkToJUNG.Label,NetworkToJUNG.Label[]> makeEdge = new Func3<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,NetworkToJUNG.Label,NetworkToJUNG.Label,NetworkToJUNG.Label[]>()
        {

            public NetworkToJUNG.Label[] execute(Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]> tree, NetworkToJUNG.Label label, NetworkToJUNG.Label label1) {

               return new NetworkToJUNG.Label[] { label, label1 };
            }
        };


        Proc2<Graph<NetworkToJUNG.Label,NetworkToJUNG.Label[]>,NetworkToJUNG.Label[]> addEdge = new Proc2<Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]>, NetworkToJUNG.Label[]>() {
            public void execute(Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]> labelGraph, NetworkToJUNG.Label[] labels) {
                labelGraph.addEdge(labels, labels[0],labels[1]);
            }
        };

         Func1<Graph<NetworkToJUNG.Label, Object>, Long> getNodeCount = new Func1<Graph<NetworkToJUNG.Label, Object>, Long>() {
            public Long execute(Graph<NetworkToJUNG.Label, Object> tree) {
                return new Long(tree.getVertexCount());
            }
        };

        Func2<Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]>, NetworkToJUNG.Label,Collection<NetworkToJUNG.Label[]>> getEdgesReachaleFromNode = new Func2<Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]>,NetworkToJUNG.Label, Collection<NetworkToJUNG.Label[]>>() {
            public Collection<NetworkToJUNG.Label[]> execute(Graph<NetworkToJUNG.Label, NetworkToJUNG.Label[]> tree, NetworkToJUNG.Label label) {
                return  new HashSet<NetworkToJUNG.Label[]>(
                        new KNeighborhoodFilter<NetworkToJUNG.Label, NetworkToJUNG.Label[]>(label, tree.getVertexCount(), KNeighborhoodFilter.EdgeType.IN_OUT).transform(tree).getEdges()
                );
            }
        };


        Proc2<Graph<NetworkToJUNG.Label, Object>,NetworkToJUNG.Label> removeNode = new Proc2<Graph<NetworkToJUNG.Label, Object>,NetworkToJUNG.Label>() {
            public void execute(Graph<NetworkToJUNG.Label, Object> tree, NetworkToJUNG.Label label) {
                tree.removeVertex(label);
            }
        };

        Func1<Graph<NetworkToJUNG.Label, Object>,NetworkToJUNG.Label> makeNode = new Func1<Graph<NetworkToJUNG.Label, Object>,NetworkToJUNG.Label>() {
            public NetworkToJUNG.Label execute(Graph<NetworkToJUNG.Label, Object> tree) {
                NetworkToJUNG.Label lbl = new NetworkToJUNG.Label(null);
                tree.addVertex(lbl);
                return lbl;
            }
        };

        Proc2<Graph<NetworkToJUNG.Label, Object>,NetworkToJUNG.Label> addNode = new Proc2<Graph<NetworkToJUNG.Label, Object>,NetworkToJUNG.Label>()
        {

            public void execute(Graph<NetworkToJUNG.Label, Object> tree, NetworkToJUNG.Label label) {
                tree.addVertex(label);
            }
        };

        TBRSearchNeighborsExpansiveInPlace tbr = new TBRSearchNeighborsExpansiveInPlace(isTreeRooted, isDestinationNode, getNodes, getEdges, getIncidentEdges, getNodesOfEdge, removeEdge, addNode, removeNode, addEdge, makeEdge, makeNode, getEdgesReachaleFromNode);

        tbr.search(inputTree, scoreTree, isScoreBetter, numIters);


        String treeNewick = JUNGToRN.toRichNewick(inputTree, isTreeRooted.execute(inputTree));
        this.richNewickGenerated(treeNewick);
        return ("\n" + treeNewick);
    }
}


