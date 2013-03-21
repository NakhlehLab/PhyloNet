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

package edu.rice.cs.bioinfo.programs.phylonet.algos.network;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import org.junit.Test;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/28/12
 * Time: 6:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneTreeProbabilityYFTest
{
    @Test
    public void test() throws Exception
    {
        BniNetwork<Integer> network = new BniNetwork();
        network.createRoot("R");

        NetNode<Integer> r = network.getRoot();
        NetNode<Integer> j = new BniNetNode<Integer>("J", null);
        NetNode<Integer> k = new BniNetNode<Integer>("K", null);
        NetNode<Integer> a = new BniNetNode<Integer>("A", null);
        NetNode<Integer> x = new BniNetNode<Integer>("X", null);
        NetNode<Integer> d = new BniNetNode<Integer>("D", null);
        NetNode<Integer> l = new BniNetNode<Integer>("L", null);
        NetNode<Integer> b = new BniNetNode<Integer>("B", null);
        NetNode<Integer> c = new BniNetNode<Integer>("C", null);

        r.adoptChild(j, 1);
        r.adoptChild(k, 1);

        k.adoptChild(d, 2);
        j.adoptChild(a, 2);
        k.adoptChild(x, 0);
        j.adoptChild(x, 0);
        x.adoptChild(l, 1);
        l.adoptChild(b, 1);
        l.adoptChild(c, 1);

        x.setParentProbability(j, .3);
        x.setParentProbability(k,.7);

        Tree gt = new STITree<Object>("(C,((B,D),A)J)K)R;");

      //  List<Double> result = new GeneTreeProbabilityYF().calculateGTDistribution(network, Arrays.asList(gt), null);



    }
}
