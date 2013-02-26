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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/19/11
 * Time: 3:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class ExtractNodeLabels implements NetworkAlgo<Collection<String>, Object, RuntimeException> {

    public ExtractNodeLabels()
    {

    }

    public Collection<String> forNetworkEmpty(NetworkEmpty network, Object input) throws RuntimeException {
        return new ArrayList<String>();
    }

    public Collection<String> forNetworkNonEmpty(NetworkNonEmpty network, Object input) throws RuntimeException {

        Collection<String> labels = new LinkedList<String>();
        collectHelp(network.PrincipleInfo, network.PrincipleDescendants, labels);

        return labels;
    }

    private void collectHelp(NetworkInfo node, DescendantList children, final Collection<String> labels) {

        node.NodeLabel.execute(new NodeLabelAlgo<Object, Object, RuntimeException>() {
            public Object forNodeLabelNonEmpty(NodeLabelNonEmpty node, Object input) throws RuntimeException {
                labels.add(node.Label.Content);
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public Object forNodeLabelEmpty(NodeLabelEmpty node, Object input) throws RuntimeException {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null);

        for(Subtree tree : children.Subtrees)
        {
            collectHelp(tree.NetworkInfo, tree.Descendants, labels);
        }

    }
}
