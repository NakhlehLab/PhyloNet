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

         blocks = Parser.parse(new ByteArrayInputStream(("#NEXUS\n" +
            "BEGIN PHYLONET;\n" +
            "\n" +
            "<z:a1,a2,a; y:b1,b2,b; c:c; d:d; e:e>;" +
            "\n" +
            "END;").getBytes()));


    }

    @Test
    public void testCommandAssigment() throws  Exception
    {
        Parser.parse(new ByteArrayInputStream(("#NEXUS\n" +
            "\n" +
            "BEGIN PHYLONET;\n" +
            "\n" +
            "\t a = SymmetricDifference network1 tree2 \"C:\\output\\symDiffResult1.txt\";\n" +
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
