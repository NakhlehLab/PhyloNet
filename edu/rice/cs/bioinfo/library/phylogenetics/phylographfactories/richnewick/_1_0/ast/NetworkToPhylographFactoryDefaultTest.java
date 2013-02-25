package edu.rice.cs.bioinfo.library.phylogenetics.phylographfactories.richnewick._1_0.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloGraph;
import edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast.RNNode;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import junit.framework.Assert;
import org.junit.Test;

import java.io.ByteArrayInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/1/12
 * Time: 2:38 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkToPhylographFactoryDefaultTest
{
    @Test
    public void testMake1()
    {
        String richNewick = "(J,K)R;";
        RichNewickReaderAST<Void> rnReader =  new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
        Networks networks = rnReader.readAnyErrorToRuntimeException(new ByteArrayInputStream(richNewick.getBytes()));

        NetworkNonEmpty networkNonEmpty = networks.Networks.iterator().next();
        PhyloGraph<RNNode> graph = new NetworkToPhyloGraphFactoryDefault().make(networkNonEmpty);

        Assert.assertEquals(3, IterableHelp.countInt(graph.getNodes()));

    }

    @Test
    public void testMake2()
    {
        String richNewick = "((A:5, (M:7)X#1:3::.7)J:1,(B:6, X#1:4::.3)K:2)R;";
        RichNewickReaderAST<Void> rnReader =  new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
        Networks networks = rnReader.readAnyErrorToRuntimeException(new ByteArrayInputStream(richNewick.getBytes()));

        NetworkNonEmpty networkNonEmpty = networks.Networks.iterator().next();
        PhyloGraph<RNNode> graph = new NetworkToPhyloGraphFactoryDefault().make(networkNonEmpty);

        Assert.assertEquals(7, IterableHelp.countInt(graph.getNodes()));

        Assert.assertEquals(1.0, graph.getEdge(new RNNode("R"), new RNNode("J")).getBranchLength());
        Assert.assertEquals(2.0, graph.getEdge(new RNNode("R"), new RNNode("K")).getBranchLength());
        Assert.assertEquals(3.0, graph.getEdge(new RNNode("J"), new RNNode("X")).getBranchLength());
        Assert.assertEquals(4.0, graph.getEdge(new RNNode("K"), new RNNode("X")).getBranchLength());
        Assert.assertEquals(5.0, graph.getEdge(new RNNode("J"), new RNNode("A")).getBranchLength());
        Assert.assertEquals(6.0, graph.getEdge(new RNNode("K"), new RNNode("B")).getBranchLength());
        Assert.assertEquals(7.0, graph.getEdge(new RNNode("X"), new RNNode("M")).getBranchLength());

        Assert.assertEquals(.7, graph.getEdge(new RNNode("J"), new RNNode("X")).getProbability());
        Assert.assertEquals(.3, graph.getEdge(new RNNode("K"), new RNNode("X")).getProbability());

    }
}
