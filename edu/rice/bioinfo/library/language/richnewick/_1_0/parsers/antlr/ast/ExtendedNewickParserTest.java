package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;

import edu.rice.bioinfo.library.language.richnewick._1_0.ast.*;
import edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.ExtendedNewickLexer;
import edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.ExtendedNewickParser;
import edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickParser;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.junit.Assert;
import org.junit.Test;

import java.io.*;
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

        Network network;



        // test the string "R;"
        network = RichNewickParser.parse(testNetworkHelp("R;"));
        AssertNodeLabelOnly("R", network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R_A;"
        // underscores in unquoted labels are replaced with space by specification
        network = RichNewickParser.parse(testNetworkHelp("R_A;"));
        AssertNodeLabelOnly("R A", network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "'R_A';"
        network = RichNewickParser.parse(testNetworkHelp("'R_A';"));
        AssertNodeLabelOnly("R_A", network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "9A9;"
        network = RichNewickParser.parse(testNetworkHelp("9A9;"));
        AssertNodeLabelOnly("9A9", network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "'Dogs'' tails wag';"
        network = RichNewickParser.parse(testNetworkHelp("'Dogs'' tails wag';"));
        AssertNodeLabelOnly("Dogs' tails wag", network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());



         // test the string "(A,B)R;"
        network = RichNewickParser.parse(testNetworkHelp("(A,B)R;"));
        AssertNodeLabelOnly("R", network.PrincipleInfo);
        Iterator<Subtree> subTrees = network.PrincipleDescendants.Subtrees.iterator();

        Subtree firstChild = subTrees.next();
        AssertNodeLabelOnly("A",firstChild.NetworkInfo);
        Assert.assertEquals(false, firstChild.Descendants.Subtrees.iterator().hasNext());

        Subtree secondChild = subTrees.next();
        AssertNodeLabelOnly("B", secondChild.NetworkInfo);
        Assert.assertEquals(false, secondChild.Descendants.Subtrees.iterator().hasNext());

          // test the string "(A,(B1,B2)B)R;"
        network = RichNewickParser.parse(testNetworkHelp("(A,(B1,B2)B)R;"));
        AssertNodeLabelOnly("R", network.PrincipleInfo);
        Iterator<Subtree> rSubTrees = network.PrincipleDescendants.Subtrees.iterator();

        Subtree rFirstChild = rSubTrees.next();
        AssertNodeLabelOnly("A",rFirstChild.NetworkInfo);
        Assert.assertEquals(false, rFirstChild.Descendants.Subtrees.iterator().hasNext());

        Subtree rSecondChild = rSubTrees.next();
        AssertNodeLabelOnly("B", rSecondChild.NetworkInfo);
        Assert.assertEquals(true, rSecondChild.Descendants.Subtrees.iterator().hasNext());

        Iterator<Subtree> bSubTrees = rSecondChild.Descendants.Subtrees.iterator();
        Subtree bFirstChild = bSubTrees.next();
        AssertNodeLabelOnly("B1",bFirstChild.NetworkInfo);
        Assert.assertEquals(false, bFirstChild.Descendants.Subtrees.iterator().hasNext());

         Subtree bSecondChild = bSubTrees.next();
        AssertNodeLabelOnly("B2",bSecondChild.NetworkInfo);
        Assert.assertEquals(false, bSecondChild.Descendants.Subtrees.iterator().hasNext());


         // test the string "R:1;"
        network = RichNewickParser.parse(testNetworkHelp("R:1;"));
        AssertNodeDetails("R", "1", null, null, null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R:1:2;"
        network = RichNewickParser.parse(testNetworkHelp("R:1:2;"));
       AssertNodeDetails("R", "1", "2", null, null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

         // test the string "R:1:2:3"
        network = RichNewickParser.parse(testNetworkHelp("R:1:2:3;"));
        AssertNodeDetails("R", "1", "2", "3", null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

         // test the string "R::2:3"
        network = RichNewickParser.parse(testNetworkHelp("R::2:3;"));
       AssertNodeDetails("R", null, "2", "3", null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R:::3"
        network = RichNewickParser.parse(testNetworkHelp("R:::3;"));
        AssertNodeDetails("R", null, null, "3", null, null, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R#0;"
        network = RichNewickParser.parse(testNetworkHelp("R#0;"));
        Assert.assertEquals("0", ((HybridNodeQualifierNonEmpty)network.PrincipleInfo.HybridNodeQualifier).HybridNodeIndex.Content);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R#H0;"
        network = RichNewickParser.parse(testNetworkHelp("R#H0;"));
        Assert.assertEquals("R", ((NodeLabelNonEmpty)network.PrincipleInfo.NodeLabel).Label.Content);
        Assert.assertEquals("0", ((HybridNodeQualifierWithType)network.PrincipleInfo.HybridNodeQualifier).HybridNodeIndex.Content);
        Assert.assertEquals("H", ((HybridNodeQualifierWithType)network.PrincipleInfo.HybridNodeQualifier).HybridNodeType.Content);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "R#H0:1:2:3;"
        network = RichNewickParser.parse(testNetworkHelp("R#H0:1:2:3;"));
        Assert.assertEquals("R", ((NodeLabelNonEmpty)network.PrincipleInfo.NodeLabel).Label.Content);
        Assert.assertEquals("0", ((HybridNodeQualifierWithType)network.PrincipleInfo.HybridNodeQualifier).HybridNodeIndex.Content);
        Assert.assertEquals("H", ((HybridNodeQualifierWithType)network.PrincipleInfo.HybridNodeQualifier).HybridNodeType.Content);
        Assert.assertEquals("1", ((BranchLengthNonEmpty)network.PrincipleInfo.BranchLength).Length.Content);
        Assert.assertEquals("2", ((BootstrapNonEmpty)network.PrincipleInfo.Bootstrap).BootstrapValue.Content);
        Assert.assertEquals("3", ((ProbabilityNonEmpty)network.PrincipleInfo.Probability).ProbabilityValue.Content);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "'This is the only node.';"
        network = RichNewickParser.parse(testNetworkHelp("'This is the only node.';"));
        AssertNodeLabelOnly("This is the only node.", network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string ";"
        network = RichNewickParser.parse(testNetworkHelp(";"));
        Assert.assertEquals(NodeLabelEmpty.Singleton, network.PrincipleInfo.NodeLabel);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        // test the string "(,);"
        network = RichNewickParser.parse(testNetworkHelp("(,);"));
        Assert.assertEquals(NodeLabelEmpty.Singleton, network.PrincipleInfo.NodeLabel);

        subTrees = network.PrincipleDescendants.Subtrees.iterator();

        firstChild = subTrees.next();
        Assert.assertEquals(NodeLabelEmpty.Singleton,firstChild.NetworkInfo.NodeLabel);
        Assert.assertEquals(false, firstChild.Descendants.Subtrees.iterator().hasNext());

        secondChild = subTrees.next();
        Assert.assertEquals(NodeLabelEmpty.Singleton, secondChild.NetworkInfo.NodeLabel);
        Assert.assertEquals(false, secondChild.Descendants.Subtrees.iterator().hasNext());

        // test the string  "((1, ((2, (3, (4)Y#H1)g)e, (((Y#H1, 5)h, 6)f)X#H2)c)a, ((X#H2, 7)d, 8)b)r;"
        network = RichNewickParser.parse(testNetworkHelp("((1, ((2, (3, (4)Y#H1)g)e, (((Y#H1, 5)h, 6)f)X#H2)c)a, ((X#H2, 7)d, 8)b)r;"));

        AssertNodeLabelOnly("r", network.PrincipleInfo);
        subTrees = network.PrincipleDescendants.Subtrees.iterator();
        Subtree rSub1 = subTrees.next();
        Subtree rSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("a", rSub1.NetworkInfo);
        subTrees = rSub1.Descendants.Subtrees.iterator();
        Subtree aSub1 = subTrees.next();
        Subtree aSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("b", rSub2.NetworkInfo);
        subTrees = rSub2.Descendants.Subtrees.iterator();
        Subtree bSub1 = subTrees.next();
        Subtree bSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("c", aSub2.NetworkInfo);
        subTrees = aSub2.Descendants.Subtrees.iterator();
        Subtree cSub1 = subTrees.next();
        Subtree cSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("d", bSub1.NetworkInfo);
        subTrees = bSub1.Descendants.Subtrees.iterator();
        Subtree dSub1 = subTrees.next();
        Subtree dSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("e", cSub1.NetworkInfo);
        subTrees = cSub1.Descendants.Subtrees.iterator();
        Subtree eSub1 = subTrees.next();
        Subtree eSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeDetails("X", null, null, null, "H", "2", dSub1.NetworkInfo);
        subTrees = dSub1.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeDetails("X", null, null, null, "H", "2", dSub1.NetworkInfo);
        subTrees = cSub2.Descendants.Subtrees.iterator();
        Subtree xSub1FromC = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("f", xSub1FromC.NetworkInfo);
        subTrees = xSub1FromC.Descendants.Subtrees.iterator();
        Subtree fSub1 = subTrees.next();
        Subtree fSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("g", eSub2.NetworkInfo);
        subTrees = eSub2.Descendants.Subtrees.iterator();
        Subtree gSub1 = subTrees.next();
        Subtree gSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeDetails("Y", null, null, null, "H", "1", gSub2.NetworkInfo);
        subTrees = gSub2.Descendants.Subtrees.iterator();
        Subtree ySub1FromG = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("h", fSub1.NetworkInfo);
        subTrees = fSub1.Descendants.Subtrees.iterator();
        Subtree hSub1 = subTrees.next();
        Subtree hSub2 = subTrees.next();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeDetails("Y", null, null, null, "H", "1", hSub1.NetworkInfo);
        subTrees = hSub1.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("1", aSub1.NetworkInfo);
        subTrees = aSub1.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("2", eSub1.NetworkInfo);
        subTrees = eSub1.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("3", gSub1.NetworkInfo);
        subTrees = gSub1.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("4", ySub1FromG.NetworkInfo);
        subTrees = ySub1FromG.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("5", hSub2.NetworkInfo);
        subTrees = hSub2.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("6", fSub2.NetworkInfo);
        subTrees = fSub2.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("7", dSub2.NetworkInfo);
        subTrees = dSub2.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());

        AssertNodeLabelOnly("8", bSub2.NetworkInfo);
        subTrees = bSub2.Descendants.Subtrees.iterator();
        Assert.assertEquals(false, subTrees.hasNext());


    }

    private void AssertNodeLabelOnly(String expectedNodeLabel, NetworkInfo info)
    {
        Assert.assertEquals(expectedNodeLabel, ((NodeLabelNonEmpty)info.NodeLabel).Label.Content);
        Assert.assertEquals(HybridNodeQualifierEmpty.Singleton, info.HybridNodeQualifier);
        Assert.assertEquals(BranchLengthEmpty.Singleton, info.BranchLength);
        Assert.assertEquals(BootstrapEmpty.Singleton, info.Bootstrap);
        Assert.assertEquals(ProbabilityEmpty.Singleton, info.Probability);
    }

    private void AssertNodeDetails(String expectedNodeLabel, String branchLength, String bootstrap, String probability, String hybridType, String hybridIndex, NetworkInfo info)
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

        if(bootstrap != null)
        {
            Assert.assertEquals(bootstrap, ((BootstrapNonEmpty)info.Bootstrap).BootstrapValue.Content);
        }
        else
        {
            Assert.assertEquals(BootstrapEmpty.Singleton, info.Bootstrap);
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
