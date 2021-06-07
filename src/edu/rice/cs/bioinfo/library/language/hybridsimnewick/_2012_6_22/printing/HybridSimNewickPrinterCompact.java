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

package edu.rice.cs.bioinfo.library.language.hybridsimnewick._2012_6_22.printing;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.io.StringWriter;
import java.util.HashSet;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/23/12
 * Time: 2:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class HybridSimNewickPrinterCompact<N> implements HybridSimNewickPrinter<N>
{

    private final Func2<N,N,String> _getBranchLength;

    private final Func2<N,N,String> _getProbability;

    private final Func1<N, String> _getHybridIndex;

    private final Func1<N, String> _getLabel;


    public HybridSimNewickPrinterCompact(Func1<N, String> getLabel, Func2<N,N,String> getBranchLength, Func2<N,N,String> getProbability,
                                         Func1<N, String> getHybridIndex)
    {
        _getLabel = getLabel;
        _getBranchLength = getBranchLength;
        _getProbability = getProbability;
        _getHybridIndex = getHybridIndex;

    }

    public void print(N root, Func1<N, Iterable<N>> getDestinationNodes, Func1<N, Tuple<N,N>> getHybridParents, String rootBranchLength, StringWriter writer)
    {
       StringBuffer buffer = writer.getBuffer();
       buffer.append(";");

       HashSet<String> seenHybridIndexes = new HashSet<String>();
       prependNode(getDestinationNodes, buffer, root, null, seenHybridIndexes, getHybridParents, rootBranchLength);
    }

    private void prependNode(Func1<N, Iterable<N>> getDestinationNodes, StringBuffer buffer, N node, N printParent, HashSet<String> seenHybridIndexes, Func1<N, Tuple<N, N>> getHybridParents, String rootBranchLength) {

        String hybridIndex = _getHybridIndex.execute(node);
        boolean isFirstVisitToHybrid = !seenHybridIndexes.contains(hybridIndex);
        if (hybridIndex != null && isFirstVisitToHybrid)
            seenHybridIndexes.add(hybridIndex);


        String nodeInfo = makeNodeInfoString(node, printParent, hybridIndex, isFirstVisitToHybrid, seenHybridIndexes, getHybridParents, rootBranchLength);
        buffer.insert(0, nodeInfo);

        if (hybridIndex != null && isFirstVisitToHybrid)
            return;


        List<N> destinationNodes = IterableHelp.toList(getDestinationNodes.execute(node));


        if (destinationNodes.size() > 0) {
            buffer.insert(0, ")");

            for (int i = 0; i < destinationNodes.size(); i++) {
                prependNode(getDestinationNodes, buffer, destinationNodes.get(i), node, seenHybridIndexes, getHybridParents, rootBranchLength);

                if (i < destinationNodes.size() - 1) {
                    buffer.insert(0, ",");
                }

            }
            buffer.insert(0, "(");
        }


    }

    private String makeNodeInfoString(N node, N printParent, String hybridIndex, boolean isFirstVisitToHybrid, HashSet<String> seenHybridIndexes, Func1<N, Tuple<N,N>> getHybridParents, String rootBranchLength) {
        String result = _getLabel.execute(node);

        if(hybridIndex != null)
        {
            if(!isFirstVisitToHybrid)
            {
                result+="#" + _getProbability.execute(printParent, node);
            }
            else
            {
                N otherHybridParent = (N) getHybridParents.execute(node).other(printParent);
                result+="#" + _getProbability.execute(otherHybridParent, node);
            }
        }

        String branchLength = printParent != null ? _getBranchLength.execute(printParent, node) : rootBranchLength;
        result += ":" + branchLength;

        return result;
    }
}

