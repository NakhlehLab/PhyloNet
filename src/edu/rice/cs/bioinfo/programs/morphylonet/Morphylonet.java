package edu.rice.cs.bioinfo.programs.morphylonet;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.*;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.Parser;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Iterator;

/**
 * Created by Matt on 5/29/2015.
 */
public class Morphylonet
{
    public static void main(String[] args)
    {

        Blocks nexusContents;
        try
        {

//            #NEXUS
//
//            Begin MORPHDATA;
//               dimensions ntax=3
//                    format symbols="1 2 3" missing=?;
//                    matrix
//            Human   1 2 3
//            Chimp   2 2 2
//            Gorilla 3 2 1;
//            END;
//
//            BEGIN PHYLONET;
//            	 SomeCommand;
//            END;
//
            FileInputStream nexusFileStream = new FileInputStream("C:\\temp\\example.nex");
            nexusContents = Parser.parse(nexusFileStream);
            Iterator<Block> blocks = nexusContents.Contents.iterator();
            MorphDataBlockBody morphData = (MorphDataBlockBody) blocks.next();
            PhyloNetBlockBody phylonet = (PhyloNetBlockBody)blocks.next();


            for(Identifier ident : morphData.MatrixContent)
            {
                System.out.print(ident.Content + " ");
            }

            System.out.println("");

            System.out.println("ntax = " + morphData.NTax.Content);
            System.out.println("symbols = " + morphData.Symbols);


        }
        catch (IOException e)
        {
            System.err.println(e.getMessage());
            System.exit(-2);
        }
        catch (CoordinateParseErrorsException e)
        {
            for(CoordinateParseError error : e.Errors)
            {
                System.err.println("Error at " + error.getLineNumber() + " " + error.getColumnNumber() + ":" + error.getMessage());
            }
            System.exit(-3);
        }

        return;


    }


}
