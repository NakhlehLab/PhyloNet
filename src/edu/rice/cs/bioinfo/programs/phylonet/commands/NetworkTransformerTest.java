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
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import junit.framework.Assert;
import org.junit.Test;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/19/11
 * Time: 5:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkTransformerTest
{
    @Test
    public void testToENewick()
    {
        NetworkInfo r = new NetworkInfo(
                new NodeLabelNonEmpty(new Text("R", -1, -1, false)),
                HybridNodeQualifierEmpty.Singleton, BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

        NetworkInfo a = new NetworkInfo(
                new NodeLabelNonEmpty(new Text("A", -1, -1, false)),
                HybridNodeQualifierEmpty.Singleton, BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

        NetworkInfo b = new NetworkInfo(
                new NodeLabelNonEmpty(new Text("B", -1, -1, false)),
                HybridNodeQualifierEmpty.Singleton, BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

         NetworkInfo x = new NetworkInfo(
                new NodeLabelNonEmpty(new Text("X", -1, -1, false)),
                HybridNodeQualifierEmpty.Singleton, BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);


         NetworkInfo z = new NetworkInfo(
                new NodeLabelNonEmpty(new Text("Z", -1, -1, false)),
                new HybridNodeQualifierNonEmpty(new Text("1", -1, -1, false)), BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

         NetworkInfo h = new NetworkInfo(
                new NodeLabelNonEmpty(new Text("H", -1, -1, false)),
                new HybridNodeQualifierNonEmpty(new Text("1", -1, -1, false)), BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

        NetworkNonEmpty oneNode = new NetworkNonEmpty(RootageQualifierEmpty.Singleton,
                                                      DescendantList.EMPTY_DESCENDANT_LIST, r);

        Assert.assertEquals("R;", NetworkTransformer.toENewick(oneNode));

        NetworkNonEmpty threeNode = new NetworkNonEmpty(RootageQualifierEmpty.Singleton, new DescendantList(a,b), r);
        Assert.assertEquals("(A,B)R;", NetworkTransformer.toENewick(threeNode));

        LinkedList<Subtree> zChildren = new LinkedList<Subtree>();
        zChildren.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, x));

        LinkedList<Subtree> hChildren = new LinkedList<Subtree>();
        hChildren.add(new Subtree(new DescendantList(zChildren), z));

        LinkedList<Subtree> aChildren = new LinkedList<Subtree>();
        aChildren.add(new Subtree(new DescendantList(hChildren), h));

        LinkedList<Subtree> bChildren = new LinkedList<Subtree>();
        bChildren.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, h));

        LinkedList<Subtree> rChildren = new LinkedList<Subtree>();
        rChildren.add(new Subtree(new DescendantList(aChildren), a));
        rChildren.add(new Subtree(new DescendantList(bChildren), b));

        NetworkNonEmpty diamond = new NetworkNonEmpty(RootageQualifierEmpty.Singleton, new DescendantList(rChildren), r);
        Assert.assertEquals("N = ((H)A,(H)B)R;\nH = (X);", NetworkTransformer.toENewick(diamond));






    }
}
