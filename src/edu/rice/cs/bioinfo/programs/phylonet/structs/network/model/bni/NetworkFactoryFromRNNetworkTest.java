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

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/15/12
 * Time: 6:00 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkFactoryFromRNNetworkTest
{
    @Test
    public void testMake()
    {

        NetworkInfo r = new NetworkInfo(new NodeLabelNonEmpty(new Text("R", -1, -1, false)), HybridNodeQualifierEmpty.Singleton, BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);
        NetworkInfo a = new NetworkInfo(new NodeLabelNonEmpty(new Text("A", -1, -1, false)), HybridNodeQualifierEmpty.Singleton, new BranchLengthNonEmpty(new Text("2", -1, -1, false)),
                                                                                             new SupportNonEmpty(new Text(".3", -1, -1, false)), ProbabilityEmpty.Singleton);
        NetworkInfo b = new NetworkInfo(new NodeLabelNonEmpty(new Text("B", -1, -1, false)), HybridNodeQualifierEmpty.Singleton, new BranchLengthNonEmpty(new Text("3", -1, -1, false)),
                                                                                             new SupportNonEmpty(new Text(".4", -1, -1, false)), ProbabilityEmpty.Singleton);

        NetworkInfo h1 = new NetworkInfo(new NodeLabelNonEmpty(new Text("H", -1, -1, false)), new HybridNodeQualifierNonEmpty(new Text("1", -1, -1, false)), new BranchLengthNonEmpty(new Text("6", -1, -1, false)),
                                                                                             new SupportNonEmpty(new Text(".7", -1, -1, false)), new ProbabilityNonEmpty(new Text(".2", -1, -1, false)));

        NetworkInfo h2 = new NetworkInfo(new NodeLabelNonEmpty(new Text("H", -1, -1, false)), new HybridNodeQualifierNonEmpty(new Text("1", -1, -1, false)), new BranchLengthNonEmpty(new Text("8", -1, -1, false)),

                                                                                              new SupportNonEmpty(new Text(".3", -1, -1, false)), new ProbabilityNonEmpty(new Text(".8", -1, -1, false)));

        Subtree aRightSub = new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, h1);
        List<Subtree> aDesc = new ArrayList();
        aDesc.add(aRightSub);

        Subtree bLeftSub = new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, h2);
        List<Subtree> bDesc = new ArrayList();
        bDesc.add(bLeftSub);


        Subtree rLeftSub = new Subtree(new DescendantList(aDesc), a);
        Subtree rRightSub = new Subtree(new DescendantList(bDesc), b);
        List<Subtree> rDesc = new ArrayList<Subtree>();
        rDesc.add(rLeftSub);
        rDesc.add(rRightSub);




        DescendantList rDescList = new DescendantList(rDesc);


        Network<Object> network = new NetworkFactoryFromRNNetwork().makeNetwork(new NetworkNonEmpty(new RootageQualifierNonEmpty(RootageQualifier.ROOTED), rDescList, r));

        NetNode<Object> netRoot = network.getRoot();

        Assert.assertEquals("R", netRoot.getName());
        Iterator<NetNode<Object>> rChildren = netRoot.getChildren().iterator();
        NetNode<Object> netA = rChildren.next();
        NetNode<Object> netB = rChildren.next();
        Assert.assertEquals("A", netA.getName());
        Assert.assertTrue(2.0 == netA.getParentDistance(netRoot));
        Assert.assertTrue(.3 == netA.getParentSupport(netRoot));
        Assert.assertTrue(1.0 == netA.getParentProbability(netRoot));
        Assert.assertEquals("B", netB.getName());
        Assert.assertTrue(3.0 == netB.getParentDistance(netRoot));
        Assert.assertTrue(.4 == netB.getParentSupport(netRoot));
        Assert.assertTrue(1.0 == netB.getParentProbability(netRoot));

        Assert.assertSame(netA.getChildren().iterator().next(), netB.getChildren().iterator().next());
        NetNode<Object> netH = netA.getChildren().iterator().next();
        Assert.assertTrue(6.0 == netH.getParentDistance(netA));
        Assert.assertTrue(.7 == netH.getParentSupport(netA));
        Assert.assertTrue(.2 == netH.getParentProbability(netA));
        Assert.assertTrue(8.0 == netH.getParentDistance(netB));
        Assert.assertTrue(.3 == netH.getParentSupport(netB));
        Assert.assertTrue(.8 == netH.getParentProbability(netB));


    }
}
