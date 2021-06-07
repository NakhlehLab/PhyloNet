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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.junit.Assert;
import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 4:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class ExtendedNewickParserTest {
    @Test
    public void testNetwork() throws Exception {

        NetworkNonEmpty network;



        // test the string "R;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R;")));
        AssertNodeLabelOnly("R", 1, 0, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R;G;"
        Networks networks = new ANTLRRichNewickParser().parse(testNetworkHelp("R;G;"));
        Iterator<NetworkNonEmpty> networkElements = networks.Networks.iterator();
        AssertNodeLabelOnly("R", 1, 0, (networkElements.next()).PrincipleInfo);
        AssertNodeLabelOnly("G", 1, 2, (networkElements.next()).PrincipleInfo);
        Assert.assertFalse(networkElements.hasNext());
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R_A;"
        // underscores in unquoted labels are replaced with space by specification
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R_A;")));
        AssertNodeLabelOnly("R A", 1, 0, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "'R_A';"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("'R_A';")));
        AssertNodeLabelOnly("R_A", 1, 0, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "9A9;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("9A9;")));
        AssertNodeLabelOnly("9A9", 1, 0, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "'Dogs'' tails wag';"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("'Dogs'' tails wag';")));
        AssertNodeLabelOnly("Dogs' tails wag", 1, 0, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());


        // test the string "(A,B)R;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("(A,B)R;")));
        AssertNodeLabelOnly("R", 1, 5, network.PrincipleInfo);
        Iterator<Subtree> subTrees = network.PrincipleDescendants.Subtrees.iterator();

        Subtree firstChild = subTrees.next();
        AssertNodeLabelOnly("A", 1, 1, firstChild.NetworkInfo);
        Assert.assertEquals(false, firstChild.Descendants.Subtrees.iterator().hasNext());

        Subtree secondChild = subTrees.next();
        AssertNodeLabelOnly("B", 1, 3, secondChild.NetworkInfo);
        Assert.assertEquals(false, secondChild.Descendants.Subtrees.iterator().hasNext());

          // test the string "(A,(B1,B2)B)R;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("(A,(B1,B2)B)R;")));
        AssertNodeLabelOnly("R", 1, 12, network.PrincipleInfo);
        Iterator<Subtree> rSubTrees = network.PrincipleDescendants.Subtrees.iterator();

        Subtree rFirstChild = rSubTrees.next();
        AssertNodeLabelOnly("A", 1, 1, rFirstChild.NetworkInfo);
        Assert.assertEquals(false, rFirstChild.Descendants.Subtrees.iterator().hasNext());

        Subtree rSecondChild = rSubTrees.next();
        AssertNodeLabelOnly("B", 1, 10, rSecondChild.NetworkInfo);
        Assert.assertEquals(true, rSecondChild.Descendants.Subtrees.iterator().hasNext());

        Iterator<Subtree> bSubTrees = rSecondChild.Descendants.Subtrees.iterator();
        Subtree bFirstChild = bSubTrees.next();
        AssertNodeLabelOnly("B1", 1, 4, bFirstChild.NetworkInfo);
        Assert.assertEquals(false, bFirstChild.Descendants.Subtrees.iterator().hasNext());

         Subtree bSecondChild = bSubTrees.next();
        AssertNodeLabelOnly("B2", 1, 7, bSecondChild.NetworkInfo);
        Assert.assertEquals(false, bSecondChild.Descendants.Subtrees.iterator().hasNext());


         // test the string "R:1;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R:1;")));
        AssertNodeDetails("R", "1", null, null, null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

          // test the string "R::2;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R::2;")));
        AssertNodeDetails("R", null, "2", null, null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

           // test the string "R:::3"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R:::3;")));
        AssertNodeDetails("R", null, null, "3", null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R:1:2;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R:1:2;")));
        AssertNodeDetails("R", "1", "2", null, null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

         // test the string "R::2:3;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R::2:3;")));
        AssertNodeDetails("R", null, "2", "3", null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R:1::3;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R:1::3;")));
        AssertNodeDetails("R", "1", null, "3", null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

         // test the string "R:1:2:3;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R:1:2:3;")));
        AssertNodeDetails("R", "1", "2", "3", null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());


        // test the string "R#0;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R#0;")));
        Assert.assertEquals("0", ((HybridNodeQualifierNonEmpty)network.PrincipleInfo.HybridNodeQualifier).HybridNodeIndex.Content);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R#H0;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R#H0;")));
        Assert.assertEquals("R", ((NodeLabelNonEmpty)network.PrincipleInfo.NodeLabel).Label.Content);
        Assert.assertEquals("0", ((HybridNodeQualifierWithType)network.PrincipleInfo.HybridNodeQualifier).HybridNodeIndex.Content);
        Assert.assertEquals("H", ((HybridNodeQualifierWithType)network.PrincipleInfo.HybridNodeQualifier).HybridNodeType.Content);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R#H0:1:2:3;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R#H0:1:2:3;")));
        Assert.assertEquals("R", ((NodeLabelNonEmpty)network.PrincipleInfo.NodeLabel).Label.Content);
        Assert.assertEquals("0", ((HybridNodeQualifierWithType)network.PrincipleInfo.HybridNodeQualifier).HybridNodeIndex.Content);
        Assert.assertEquals("H", ((HybridNodeQualifierWithType)network.PrincipleInfo.HybridNodeQualifier).HybridNodeType.Content);
        Assert.assertEquals("1", ((BranchLengthNonEmpty)network.PrincipleInfo.BranchLength).Length.Content);
        Assert.assertEquals("2", ((SupportNonEmpty)network.PrincipleInfo.Support).SupportValue.Content);
        Assert.assertEquals("3", ((ProbabilityNonEmpty)network.PrincipleInfo.Probability).ProbabilityValue.Content);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "[&R]R;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("[&R]R;")));
        AssertNodeLabelOnly("R", 1, 4, network.PrincipleInfo);
        Assert.assertEquals("[&R]", ((RootageQualifierNonEmpty)network.RootageQualifier).Qualifier);

         // test the string "R;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("R;")));
        AssertNodeLabelOnly("R", 1, 0, network.PrincipleInfo);
        Assert.assertEquals(RootageQualifierEmpty.Singleton, network.RootageQualifier);

        // test the string "'This is the only node.';"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("'This is the only node.';")));
        AssertNodeLabelOnly("This is the only node.", 1, 0, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string ";"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp(";")));
        Assert.assertEquals(NodeLabelEmpty.Singleton, network.PrincipleInfo.NodeLabel);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "(,);"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("(,);")));
        Assert.assertEquals(NodeLabelEmpty.Singleton, network.PrincipleInfo.NodeLabel);

        subTrees = network.PrincipleDescendants.Subtrees.iterator();

        firstChild = subTrees.next();
        Assert.assertEquals(NodeLabelEmpty.Singleton,firstChild.NetworkInfo.NodeLabel);
        Assert.assertEquals(false, firstChild.Descendants.Subtrees.iterator().hasNext());

        secondChild = subTrees.next();
        Assert.assertEquals(NodeLabelEmpty.Singleton, secondChild.NetworkInfo.NodeLabel);
        Assert.assertEquals(false, secondChild.Descendants.Subtrees.iterator().hasNext());

        // test the string  "((1, ((2, (3, (4)Y#H1)g)e, (((Y#H1, 5)h, 6)f)X#H2)c)a, ((X#H2, 7)d, 8)b)r;"
        network = SingleNetwork(new ANTLRRichNewickParser().parse(testNetworkHelp("((1, ((2, (3, (4)Y#H1)g)e, (((Y#H1, 5)h, 6)f)X#H2)c)a, ((X#H2, 7)d, 8)b)r;")));

        AssertNodeLabelOnly("r", 1, 72, network.PrincipleInfo);
        subTrees = network.PrincipleDescendants.Subtrees.iterator();
        Subtree rSub1 = subTrees.next();
        Subtree rSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("a", 1, 52, rSub1.NetworkInfo);
        subTrees = rSub1.Descendants.Subtrees.iterator();
        Subtree aSub1 = subTrees.next();
        Subtree aSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("b", 1, 70, rSub2.NetworkInfo);
        subTrees = rSub2.Descendants.Subtrees.iterator();
        Subtree bSub1 = subTrees.next();
        Subtree bSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("c", 1, 50, aSub2.NetworkInfo);
        subTrees = aSub2.Descendants.Subtrees.iterator();
        Subtree cSub1 = subTrees.next();
        Subtree cSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("d", 1, 65, bSub1.NetworkInfo);
        subTrees = bSub1.Descendants.Subtrees.iterator();
        Subtree dSub1 = subTrees.next();
        Subtree dSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("e", 1, 24, cSub1.NetworkInfo);
        subTrees = cSub1.Descendants.Subtrees.iterator();
        Subtree eSub1 = subTrees.next();
        Subtree eSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeDetails("X", null, null, null, "H", "2", dSub1.NetworkInfo);
        subTrees = dSub1.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeDetails("X",  null, null, null, "H", "2", dSub1.NetworkInfo);
        subTrees = cSub2.Descendants.Subtrees.iterator();
        Subtree xSub1FromC = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("f", 1, 43, xSub1FromC.NetworkInfo);
        subTrees = xSub1FromC.Descendants.Subtrees.iterator();
        Subtree fSub1 = subTrees.next();
        Subtree fSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("g", 1, 22, eSub2.NetworkInfo);
        subTrees = eSub2.Descendants.Subtrees.iterator();
        Subtree gSub1 = subTrees.next();
        Subtree gSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeDetails("Y", null, null, null, "H", "1", gSub2.NetworkInfo);
        subTrees = gSub2.Descendants.Subtrees.iterator();
        Subtree ySub1FromG = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("h", 1, 38, fSub1.NetworkInfo);
        subTrees = fSub1.Descendants.Subtrees.iterator();
        Subtree hSub1 = subTrees.next();
        Subtree hSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeDetails("Y", null, null, null, "H", "1", hSub1.NetworkInfo);
        subTrees = hSub1.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("1", 1, 2, aSub1.NetworkInfo);
        subTrees = aSub1.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("2", 1, 7, eSub1.NetworkInfo);
        subTrees = eSub1.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("3", 1, 11, gSub1.NetworkInfo);
        subTrees = gSub1.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("4", 1, 15, ySub1FromG.NetworkInfo);
        subTrees = ySub1FromG.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("5", 1, 36, hSub2.NetworkInfo);
        subTrees = hSub2.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("6", 1, 41, fSub2.NetworkInfo);
        subTrees = fSub2.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("7", 1, 63, dSub2.NetworkInfo);
        subTrees = dSub2.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("8", 1, 68, bSub2.NetworkInfo);
        subTrees = bSub2.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());


    }

    private NetworkNonEmpty SingleNetwork(Networks networks)
    {
       Iterator<NetworkNonEmpty> i = networks.Networks.iterator();

       NetworkNonEmpty first = (NetworkNonEmpty)i.next();

       Assert.assertFalse(i.hasNext());

       return first;

    }

    private void AssertNodeLabelOnly(String expectedNodeLabel, int lineNumber, int columnNumber, NetworkInfo info)
    {
        NodeLabelNonEmpty nodeLabel =    ((NodeLabelNonEmpty)info.NodeLabel);

        Assert.assertEquals(expectedNodeLabel, nodeLabel.Label.Content);
        Assert.assertEquals(lineNumber, nodeLabel.Label.LineNumberStart);
        Assert.assertEquals(columnNumber, nodeLabel.Label.ColumnNumberStart);
        Assert.assertEquals(HybridNodeQualifierEmpty.Singleton, info.HybridNodeQualifier);
        Assert.assertEquals(BranchLengthEmpty.Singleton, info.BranchLength);
        Assert.assertEquals(SupportEmpty.Singleton, info.Support);
        Assert.assertEquals(ProbabilityEmpty.Singleton, info.Probability);
    }

    private void AssertNodeDetails(String expectedNodeLabel, String branchLength, String support, String probability, String hybridType, String hybridIndex, NetworkInfo info)
    {
        if(expectedNodeLabel != null)
        {
            Assert.assertEquals(expectedNodeLabel, ((NodeLabelNonEmpty)info.NodeLabel).Label.Content);
        }
        else
        {
            Assert.assertEquals(NodeLabelEmpty.Singleton, info.NodeLabel);
        }

        if(branchLength != null)
        {
            Assert.assertEquals(branchLength, ((BranchLengthNonEmpty)info.BranchLength).Length.Content);
        }
        else
        {
            Assert.assertEquals(BranchLengthEmpty.Singleton, info.BranchLength);
        }

        if(support != null)
        {
            Assert.assertEquals(support, ((SupportNonEmpty)info.Support).SupportValue.Content);
        }
        else
        {
            Assert.assertEquals(SupportEmpty.Singleton, info.Support);
        }

        if(probability != null)
        {
            Assert.assertEquals(probability, ((ProbabilityNonEmpty)info.Probability).ProbabilityValue.Content);
        }
        else
        {
            Assert.assertEquals(ProbabilityEmpty.Singleton, info.Probability);
        }

        if(hybridType != null)
        {
            Assert.assertEquals(hybridType, ((HybridNodeQualifierWithType)info.HybridNodeQualifier).HybridNodeType.Content);
        }

        if(hybridIndex != null)
        {
            Assert.assertEquals(hybridIndex, ((HybridNodeQualifierNonEmpty)info.HybridNodeQualifier).HybridNodeIndex.Content);
        }
        else
        {
            Assert.assertEquals(HybridNodeQualifierEmpty.Singleton, info.HybridNodeQualifier);
        }


    }

    private ExtendedNewickParser testNetworkHelp(String network) throws Exception
    {
        ANTLRInputStream inStream = new ANTLRInputStream(new ByteArrayInputStream(network.getBytes()));
        ExtendedNewickLexer lexer = new ExtendedNewickLexer(inStream);
        return new ExtendedNewickParser(new CommonTokenStream(lexer));
    }

}
