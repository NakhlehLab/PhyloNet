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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.GraphBuilder;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseGraph;

import java.math.BigDecimal;
import java.math.BigInteger;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/23/12
 * Time: 5:07 PM
 * To change this template use File | Settings | File Templates.
 */
class NetworkToJUNG
{
    public static class Label
    {
        public final String Content;

        public Label(String content)
        {
            Content = content;
        }

        public String toString()
        {
            return Content == null ? "" : Content;
        }
    }


    public Graph<Label,Label[]> toTreeWithoutEdgeProperties(NetworkNonEmpty net)
    {

       final Graph<Label,Label[]> result = net.RootageQualifier.execute(new RootageQualifierAlgo<Graph<Label,Label[]>, Object, RuntimeException>() {
            public Graph<Label,Label[]> forEmptyQualifier(RootageQualifierEmpty rootageQualifierEmpty, Object o) throws RuntimeException {
                return new DirectedSparseGraph<Label, Label[]>();
            }

            public Graph<Label,Label[]> forNonEmptyQualifier(RootageQualifierNonEmpty rootageQualifierNonEmpty, Object o) throws RuntimeException {

                if(rootageQualifierNonEmpty.isRooted())
                {
                    return new DirectedSparseGraph<Label, Label[]>();
                }
                return new SparseGraph<Label, Label[]>();
            }
        }, null);

        GraphBuilder<Label> builder = new GraphBuilder<Label>() {
            public Label createNode(String label) {
               Label lbl = new Label(label);
                result.addVertex(lbl);
                return lbl;

            }

            public Label createHybridNode(String label, HybridNodeType hybridNodeType, BigInteger hybridNodeIndex) {
                throw new IllegalArgumentException("Passed network should be a tree, but found hybrid node.");
            }

            public void createDirectedEdge(Label node1, Label node2, BigDecimal bigDecimal, BigDecimal bigDecimal1, BigDecimal bigDecimal2) {
                // unrooted trees have no directed edge, so we will ignore the direction
                // we also ignore any edge properties on the tree we construct

                Label[] edge = {node1, node2};
                result.addEdge(edge, node1, node2);


            }
        };

        DAGFactory.makeDAG(net, builder);

        return result;
    }


}
