package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.SingleLinePrinter;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickReaderAST_ANTLR;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.StringReader;
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
        RichNewickReaderAST_ANTLR reader = new RichNewickReaderAST_ANTLR();

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
