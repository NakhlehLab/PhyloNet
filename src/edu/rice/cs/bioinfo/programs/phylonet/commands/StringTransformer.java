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

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.SingleLinePrinter;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/31/11
 * Time: 2:44 PM
 * To change this template use File | Settings | File Templates.
 */
class StringTransformer
{
    static String toRNewickTree(String eNewickTree) throws CoordinateParseErrorsException
    {
        RichNewickReaderAST reader =  new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);

        try
        {
           RichNewickReadResult<Networks> result = reader.read(new ByteArrayInputStream(eNewickTree.getBytes()));
           Iterator<NetworkNonEmpty> networks = result.getNetworks().Networks.iterator();
           NetworkNonEmpty tree = networks.next();

           if(networks.hasNext())
           {
               throw new IllegalArgumentException("Passed eNewick string was not a single tree.");
           }

           SingleLinePrinter printer = new SingleLinePrinter();
           printer.setSupportTransformer(new TransformSupportToBase100());
           return printer.toString(tree);

        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }
}
