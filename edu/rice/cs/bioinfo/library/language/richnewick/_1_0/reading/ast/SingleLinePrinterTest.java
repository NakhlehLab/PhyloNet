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

import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 5:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class SingleLinePrinterTest {

    @Test
    public void testToString() {
        NetworkInfo a = new NetworkInfo(new NodeLabelNonEmpty(new Text("A", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);
        NetworkInfo a1 = new NetworkInfo(new NodeLabelNonEmpty(new Text("A1", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);
        NetworkInfo a2 = new NetworkInfo(new NodeLabelNonEmpty(new Text("A2", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

        NetworkInfo b = new NetworkInfo(new NodeLabelNonEmpty(new Text("B", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);
        NetworkInfo c = new NetworkInfo(new NodeLabelNonEmpty(new Text("C", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);
        NetworkInfo r = new NetworkInfo(new NodeLabelNonEmpty(new Text("R", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

        ArrayList<Subtree> aChildren = new ArrayList<Subtree>();
        aChildren.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, a1));
        aChildren.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, a2));

        ArrayList<Subtree> networkDl = new ArrayList<Subtree>();
        networkDl.add(new Subtree(new DescendantList(aChildren), a));
        networkDl.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, b));
        networkDl.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, c));


        NetworkNonEmpty network = new NetworkNonEmpty(RootageQualifierEmpty.Singleton, new DescendantList(networkDl), r);

        Assert.assertEquals("((A1,A2)A,B,C)R;", new SingleLinePrinter().toString(network));

         NetworkInfo z = new NetworkInfo(new NodeLabelNonEmpty(new Text("Z", 1, 0, false)), new HybridNodeQualifierNonEmpty(new Text("1", -1, -1, false)),
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

        ArrayList<Subtree> a1Children = new ArrayList<Subtree>();
        a1Children.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, z));

        ArrayList<Subtree> a2Children = new ArrayList<Subtree>();
        a2Children.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, z));

        aChildren = new ArrayList<Subtree>();
        aChildren.add(new Subtree(new DescendantList(a1Children), a1));
        aChildren.add(new Subtree(new DescendantList(a2Children), a2));

        network = new NetworkNonEmpty(RootageQualifierEmpty.Singleton, new DescendantList(aChildren), a);

        Assert.assertEquals("((Z#1)A1,(Z#1)A2)A;", new SingleLinePrinter().toString(network));


    }
}
