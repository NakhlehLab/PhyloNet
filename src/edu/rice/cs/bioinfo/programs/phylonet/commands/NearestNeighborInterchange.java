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
import edu.rice.cs.bioinfo.programs.phylonet.algos.nni.NNISearcherNeighborsExpansiveInPlace;
import edu.uci.ics.jung.graph.Graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/21/12
 * Time: 2:29 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("nninterchange")
public class NearestNeighborInterchange extends HeuristicTreeSearchBase {


     public NearestNeighborInterchange(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
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


        NNISearcherNeighborsExpansiveInPlace nni = new NNISearcherNeighborsExpansiveInPlace(isTreeRooted, isDestinationNode, getNodes, getEdges,  getIncidentEdges, getNodesOfEdge, makeEdge, addEdge, removeEdge);


        nni.search(inputTree, scoreTree, isScoreBetter, numIters);

        String treeNewick = JUNGToRN.toRichNewick(inputTree, isTreeRooted.execute(inputTree));

        this.richNewickGenerated(treeNewick);
        return("\n" + treeNewick);
    }
}
