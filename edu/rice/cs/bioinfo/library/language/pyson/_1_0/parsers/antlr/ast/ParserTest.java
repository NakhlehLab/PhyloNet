package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.Blocks;
import junit.framework.Assert;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.junit.Test;

import java.io.ByteArrayInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/1/11
 * Time: 5:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParserTest {

    @Test
    public void testBlocks() throws Exception
    {
        Blocks blocks;

        blocks = Parser.parse(new ByteArrayInputStream("#NEXUS".getBytes()));
        Assert.assertFalse(blocks.Contents.iterator().hasNext());

        blocks = Parser.parse(new ByteArrayInputStream(("#NEXUS\n" +
            "\n" +
            "BEGIN TREES;\n" +
            "\n" +
            "\tTree tree1 = ((1,2),3);\n" +
            "\tTree tree2 = ((beetle,fly),spider);\n" +
            "\tTree tree3 = ((Scarabaeus,Drosophila),Aranaeus);\n" +
            "\n" +
            "END;\n" +
            "\n" +
            "BEGIN PHYLONET;\n" +
            "\n" +
            "\tSymmetricDifference network1 tree2 \"C:\\output\\symDiffResult1.txt\";\n" +
            "\tSymmetricDifference network2 tree3 \"C:\\output\\symDiffResult2.txt\";\n" +
            "\n" +
            "END;").getBytes()));


    }

    private PySONParser testBlocksHelp(String network) throws Exception
    {
        ANTLRInputStream inStream = new ANTLRInputStream(new ByteArrayInputStream(network.getBytes()));
        PySONLexer lexer = new PySONLexer(inStream);
        return new PySONParser(new CommonTokenStream(lexer));
    }
}
