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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/26/11
 * Time: 1:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphBuilderDOT implements GraphBuilder<Integer>
{
    private int _nextNodeNumber = 1;

    private  StringBuffer _graph = new StringBuffer("digraph network {");

    public Integer createNode(String label)
    {
        try
        {
            String dotLabel = label == null ? "" : label;
            _graph.append("\n\t" + _nextNodeNumber + "[label=\"" + dotLabel + "\"];");
            return new Integer(_nextNodeNumber);
        }
        finally
        {
            _nextNodeNumber++;
        }

    }

    public Integer createHybridNode(String label, HybridNodeType hybridType, BigInteger hybridNodeIndex) {

        try
        {
            String dotLabel = "(" + hybridNodeIndex + (makeHybridTypeLabel(hybridType)) + ")";

            if(label != null)
                dotLabel = label + " " + dotLabel;

            _graph.append("\n\t" + _nextNodeNumber + "[label=\"" + dotLabel + "\" shape=\"box\"];");
            return new Integer(_nextNodeNumber);
        }
        finally
        {
            _nextNodeNumber++;
        }
    }

    private String makeHybridTypeLabel(HybridNodeType hybridType) {

        if(hybridType == null || hybridType == HybridNodeType.Unspecified)
            return "";

        if(hybridType == HybridNodeType.LateralGeneTransfer)
            return "LGT";

        if(hybridType == HybridNodeType.Hybridization)
            return "H";

        if(hybridType == HybridNodeType.Recombination)
            return "R";

         throw new IllegalArgumentException("Unexpected hybrid type.");
    }

    public void createDirectedEdge(Integer tail, Integer tip, BigDecimal branchLength, BigDecimal support, BigDecimal probability) {

        String label = null;

        if(branchLength != null || support != null || probability != null)
        {
            label = "[label=\"";

            ArrayList<String> adj = new ArrayList<String>();

            if(branchLength != null)
                adj.add("bl:" + branchLength.toPlainString());

            if(support != null)
                adj.add("bs:" + support.toPlainString());

            if(probability != null)
                adj.add("p:" + probability);

            for(int i = 0; i<adj.size(); i++)
            {
                label+= adj.get(i).replace('"', '\'');

                if(i != adj.size() -1)
                    label+=" ";
            }


            label+= "\"]";
        }

        _graph.append("\n\t" + tail + "->" + tip + (label == null ? "" : label) + ";");

    }

    public String getDOT()
    {
        _graph.append("\n}");
        return _graph.toString();
    }
}
