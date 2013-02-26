/*
 * Copyright (c) 2013 Rice University.
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

package edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.csa;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.IsRooted;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkInfo;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.csa.ASTNetworkInspector;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.csa.CSAError;
import edu.rice.cs.bioinfo.library.programming.Func1;

import java.math.BigDecimal;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/21/13
 * Time: 5:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class ASTContextAnalyser
{
    public Iterable<CSAError> analyse(Network network, final BigDecimal hybridSumTollerance)
    {
        NetworkNonEmpty nonEmptyNetwork = network.execute(new NetworkAlgo<NetworkNonEmpty, RuntimeException>() {
            public NetworkNonEmpty forNetworkEmpty(NetworkEmpty network) {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public NetworkNonEmpty forNetworkNonEmpty(NetworkNonEmpty network) {
                return network;  //To change body of implemented methods use File | Settings | File Templates.
            }
        });

        if(nonEmptyNetwork != null)
        {

            final ASTNetworkInspector inspector = new ASTNetworkInspector(network);
            final Func1<Object,NetworkInfo> networkNodeToPrimarySyntaxNode = new Func1<Object, NetworkInfo>() {
                public NetworkInfo execute(Object networkNode) {

                    return inspector.getPrimarySyntaxNode(networkNode);
                }
            };

            final boolean isRooted = nonEmptyNetwork.RootageQualifier.execute(new IsRooted(), null);
            final ContextAnalyser csAnalyser = new ContextAnalyser(hybridSumTollerance);

            return nonEmptyNetwork.TreeProbability.execute(new TreeProbabilityAlgo<Iterable<CSAError>, Exception>() {
                public Iterable<CSAError> forEmpty(TreeProbabilityEmpty empty) {
                    return csAnalyser.analyse(inspector.getSyntaxNodes(), inspector, inspector.getNetworkNodes(),
                                              inspector, networkNodeToPrimarySyntaxNode, isRooted);
                }

                public Iterable<CSAError> forNonEmpty(TreeProbabilityNonEmpty nonEmpty) {
                    return csAnalyser.analyse(inspector.getSyntaxNodes(), inspector, inspector.getNetworkNodes(),
                                              inspector, networkNodeToPrimarySyntaxNode, isRooted, nonEmpty.ProbString, nonEmpty.LineNumber, nonEmpty.ColumnNumber);
                }
            });

        }

        return new ArrayList<CSAError>();
    }
}





