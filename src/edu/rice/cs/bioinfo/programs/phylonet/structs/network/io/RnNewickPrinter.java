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

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.io;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.RichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.io.StringWriter;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/1/12
 * Time: 12:53 PM
 * To change this template use File | Settings | File Templates.
 */
public class RnNewickPrinter<T> extends RichNewickPrinterCompact<NetNode<T>>
{
    private final Func1<NetNode<T>, String> _getLabel = new Func1<NetNode<T>, String>() {
            public String execute(NetNode<T> input) {
                return input.getName().replace(' ', '_');
            }
        };

    private final Func1<NetNode<T>, Iterable<NetNode<T>>> _getDestinationNodes = new Func1<NetNode<T>, Iterable<NetNode<T>>>() {
            public Iterable<NetNode<T>> execute(NetNode<T> input) {
                return input.getChildren();
            }
        };

    private final Func1<NetNode<T>,HybridNodeType> _getHybridNodeType = new Func1<NetNode<T>, HybridNodeType>() {
        public HybridNodeType execute(NetNode<T> input) {
            return HybridNodeType.Hybridization;
        }
    };

    public RnNewickPrinter()
    {
        this.setGetBranchLength(new Func2<NetNode<T>, NetNode<T>, String>() {
            public String execute(NetNode<T> parent, NetNode<T> child) {
                double parentDistance = child.getParentDistance(parent);
                if(NetNode.NO_DISTANCE == parentDistance)
                {
                    return null;
                }
                return  parentDistance + "";
            }
        });

        this.setGetProbability(new Func2<NetNode<T>, NetNode<T>, String>() {
            public String execute(NetNode<T> parent, NetNode<T> child) {

                if(child.getIndeg() < 2)
                    return null;

                double probability = child.getParentProbability(parent);
                return probability==NetNode.NO_PROBABILITY ? null : probability + "";

            }
        });

        this.setGetSupport(new Func2<NetNode<T>, NetNode<T>, String>() {
            public String execute(NetNode<T> parent, NetNode<T> child) {
                double support = child.getParentSupport(parent);
                return support==NetNode.NO_SUPPORT ? null : support + "";
            }
        });
    }

    public void print(Network<T> network, StringWriter writer)
    {

        final Map<NetNode<T>, Integer> hybridNodeToHybridIndex = new HashMap<NetNode<T>, Integer>();

        for(NetNode<T> node : network.dfs())
        {
            if(node.getIndeg() > 1)
            {
                hybridNodeToHybridIndex.put(node, hybridNodeToHybridIndex.size() + 1);
            }
        }

        Func1<NetNode<T>, String> getHybridIndex = new Func1<NetNode<T>, String>() {
            public String execute(NetNode<T> input) {
                if(hybridNodeToHybridIndex.containsKey(input))
                {
                    return hybridNodeToHybridIndex.get(input).toString();
                }
                else
                {
                    return null;
                }
            }
        };


        print(true, network.getRoot(), _getLabel, _getDestinationNodes, getHybridIndex, _getHybridNodeType, writer);



    }
}
