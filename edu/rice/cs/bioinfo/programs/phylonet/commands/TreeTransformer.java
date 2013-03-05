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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/27/11
 * Time: 1:20 PM
 * To change this template use File | Settings | File Templates.
 */
class TreeTransformer {

    static Network toNetwork(Tree tree)
    {
        if(tree.isEmpty())
        {
            return NetworkEmpty.Singleton;
        }
        else if(tree.isRooted())
        {
            Subtree network = makeSubTree(tree.getRoot());

            return new NetworkNonEmpty(new RootageQualifierNonEmpty(RootageQualifier.ROOTED), network.Descendants, network.NetworkInfo);
        }
        else
        {
            Subtree network = makeSubTree(tree.getNodes().iterator().next());
            return new NetworkNonEmpty(new RootageQualifierNonEmpty(RootageQualifier.UNROOTED), network.Descendants, network.NetworkInfo);
        }

    }

    private static Subtree makeSubTree(TNode node)
    {
        NetworkInfo nodeInfo = makeNetworkInfo(node);
        DescendantList descendantsList = makeDescendantsList(node);

        return new Subtree(descendantsList, nodeInfo);

    }

    private static NetworkInfo makeNetworkInfo(TNode node)
    {
        String nodeName = node.getName();

        NodeLabel nodeLabel = nodeName == TNode.NO_NAME ?
                                NodeLabelEmpty.Singleton :
                                new NodeLabelNonEmpty(new Text(nodeName, -1, -1, false));

        return new NetworkInfo(nodeLabel, HybridNodeQualifierEmpty.Singleton, BranchLengthEmpty.Singleton, SupportEmpty.Singleton,
                               ProbabilityEmpty.Singleton);
    }

    private static DescendantList makeDescendantsList(TNode node) {

        LinkedList<Subtree> childTrees = new LinkedList<Subtree>();

        for(TNode child : node.getChildren())
        {
            Subtree childTree = makeSubTree(child);
            childTrees.addLast(childTree);
        }

        return new DescendantList(childTrees);
    }
}
