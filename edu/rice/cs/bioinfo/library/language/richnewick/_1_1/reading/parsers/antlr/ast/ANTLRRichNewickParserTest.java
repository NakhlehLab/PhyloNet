/*
 * Copyright (c) 2013 Rice University.
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

package edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickParser;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/22/13
 * Time: 2:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class ANTLRRichNewickParserTest extends RichNewickParserOfNetworksTest<Networks> {

    @Override
    protected RichNewickParser<Networks> makeParser() {

        return new RichNewickParser<Networks>() {
            public Networks parse(final InputStream instream) throws CoordinateParseErrorsException {

                try
                {
                    return new ANTLRRichNewickParser().parse(instream);
                }
                catch(IOException e)
                {
                    throw new RuntimeException(e);
                }

            }
        };
    }

    @Test
    public void testNetworksWithTreeProb()  throws Exception
    {
        // test the string "[&W .9]R;"
        NetworkNonEmpty network = SingleNetwork(makeParser().parse(toInputStream("[&W .9]R;")));
        AssertNodeLabelOnly("R", 1, 7, network.PrincipleInfo);
        Assert.assertEquals(false, network.PrincipleDescendants.Subtrees.iterator().hasNext());

        network.TreeProbability.execute(new TreeProbabilityAlgo<Void, RuntimeException>() {

            public Void forEmpty(TreeProbabilityEmpty empty) {
                Assert.fail();
                return null;
            }

            public Void forNonEmpty(TreeProbabilityNonEmpty nonEmpty) {
               Assert.assertEquals("[&W .9]", nonEmpty.Content);
               return null;
            }
        });
    }

    private NetworkNonEmpty SingleNetwork(Networks networks)
    {
        Iterator<NetworkNonEmpty> i = networks.Networks.iterator();

        NetworkNonEmpty first = (NetworkNonEmpty)i.next();

        Assert.assertFalse(i.hasNext());

        return first;

    }
}
