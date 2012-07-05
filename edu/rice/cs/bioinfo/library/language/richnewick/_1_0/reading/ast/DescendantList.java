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
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 6:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class DescendantList implements AbstractSyntaxNode
{
    public static final DescendantList EMPTY_DESCENDANT_LIST = new DescendantList(new ArrayList<Subtree>());

    public final Iterable<Subtree> Subtrees;

    public DescendantList(Iterable<Subtree> subtrees)
    {
        Subtrees = subtrees;
    }

    public DescendantList(NetworkInfo... leafs)
    {
        LinkedList<Subtree> trees = new LinkedList<Subtree>();

        for(NetworkInfo n : leafs)
        {
            trees.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, n));
        }

        Subtrees = trees;
    }
}
