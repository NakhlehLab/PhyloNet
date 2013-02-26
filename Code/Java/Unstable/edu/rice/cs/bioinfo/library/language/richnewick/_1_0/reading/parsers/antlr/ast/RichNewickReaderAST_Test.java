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

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.GraphBuilderNoAction;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import org.junit.Assert;
import org.junit.Test;

import java.io.ByteArrayInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/1/11
 * Time: 5:34 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RichNewickReaderAST_Test {

    protected abstract RichNewickReaderAST makeReader();

    @Test
    public void testRead() throws Exception
    {
        RichNewickReaderAST reader = makeReader();

        try
        {
            reader.read(new ByteArrayInputStream("((H#1:::.6)A,(H#1:::.6)B)C;".getBytes()), GraphBuilderNoAction.Singleton);
            Assert.fail("Expected exception.");
        }
        catch(CoordinateParseErrorsException e)
        {

        }

           reader.read(
                    new ByteArrayInputStream("((1, ((2, (3, (4)Y#H1)g)e, (((Y#H1, 5)h, 6)f)X#H2)c)a, ((X#H2, 7)d, 8)b)r;".getBytes()),
                    GraphBuilderNoAction.Singleton);



        RichNewickReadResult<Networks> result = reader.read(
            new ByteArrayInputStream("((X#H1::0.5)A,(X#H1::0.5)B)R;".getBytes()),
            GraphBuilderNoAction.Singleton);
        NetworkNonEmpty rNode = result.getNetworks().Networks.iterator().next();
        Subtree sTree = rNode.PrincipleDescendants.Subtrees.iterator().next();
        Subtree xTree = sTree.Descendants.Subtrees.iterator().next();
        Assert.assertEquals("X", ((NodeLabelNonEmpty)xTree.NetworkInfo.NodeLabel).Label.Content);
        Assert.assertEquals("0.5", ((ProbabilityNonEmpty)xTree.NetworkInfo.Probability).ProbabilityValue.Content);



    }
}
