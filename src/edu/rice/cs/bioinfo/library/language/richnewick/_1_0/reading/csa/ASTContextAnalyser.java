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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.csa;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.csa.CSAError;
import edu.rice.cs.bioinfo.library.programming.Func1;

import java.math.BigDecimal;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/8/11
 * Time: 2:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class ASTContextAnalyser {

    public static Iterable<CSAError> analyse(Network network, BigDecimal hybridSumTollerance)
    {
         final ASTNetworkInspector inspector = new ASTNetworkInspector(network);
         Func1<Object,NetworkInfo> networkNodeToPrimarySyntaxNode = new Func1<Object, NetworkInfo>() {
               public NetworkInfo execute(Object networkNode) {

                    return inspector.getPrimarySyntaxNode(networkNode);
                }
               };

        boolean isRooted = network.execute(new NetworkAlgo<Boolean, Object, RuntimeException>() {
            public Boolean forNetworkEmpty(NetworkEmpty network, Object input) throws RuntimeException {
                return false;
            }

            public Boolean forNetworkNonEmpty(NetworkNonEmpty network, Object input) throws RuntimeException {

                return network.RootageQualifier.execute(new RootageQualifierAlgo<Boolean, Object, RuntimeException>() {

                    public Boolean forEmptyQualifier(RootageQualifierEmpty rootage, Object input) throws RuntimeException {
                       return true;
                    }

                    public Boolean forNonEmptyQualifier(RootageQualifierNonEmpty rootage, Object input) throws RuntimeException {

                        return rootage.isRooted();
                    }
                }, null);
            }
        }, null);

        return new ContextAnalyser(hybridSumTollerance).analyse(inspector.getSyntaxNodes(), inspector, inspector.getNetworkNodes(),
                                       inspector, networkNodeToPrimarySyntaxNode, isRooted);
    }
}
