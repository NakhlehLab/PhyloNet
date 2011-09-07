package edu.rice.bioinfo.programs.phylonet;

import edu.rice.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.bioinfo.library.language.pyson._1_0.ast.Blocks;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.SACFactoryFromAST;
import edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.*;
import edu.rice.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands.StringsAndCommands;
import edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.PySONLexer;
import edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.PySONParser;
import edu.rice.bioinfo.library.programming.Proc3;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;

import javax.swing.text.html.parser.Parser;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 1:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program {

    public static void main(String[] args) throws FileNotFoundException, IOException
    {
        if(args.length != 1)
        {
            showUsage();
            return;
        }

        File inFile = new File(args[0]);

        if(!inFile.exists())
        {
            showFileDoesNotExist(inFile);
        }

        try
        {
            Blocks blocks = parseToBlocks(new FileInputStream(inFile));
            StringsAndCommands sac = SACFactoryFromAST.make(blocks);

            ContextSensitiveAnalyser.Analyse(sac, new Proc3<String, Integer, Integer>() {

                public void execute(String message, Integer line, Integer col) {

                    System.err.println(String.format("Error at [%s,%s]: %s", line, col, message));

                }
            });
        }
        catch(CoordinateParseErrorsException e)
        {
            showParseError(e);
            return;
        }




    }

    private static void showParseError(CoordinateParseErrorsException e) {
    }

    private static Blocks parseToBlocks(FileInputStream fileInputStream) throws IOException, CoordinateParseErrorsException {

        return edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.Parser.parse(fileInputStream);

    }

    private static void showFileDoesNotExist(File inFile) {
        //To change body of created methods use File | Settings | File Templates.
    }

    private static void showUsage() {
        //To change body of created methods use File | Settings | File Templates.
    }
}
