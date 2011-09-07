package edu.rice.bioinfo.programs.pysonparse;

import edu.rice.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.Parser;
import org.antlr.runtime.RecognitionException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/1/11
 * Time: 5:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program {

    public static void main(String[] args) throws Exception
    {
        if(args.length < 1)
        {
            System.err.println("Usage: java pysonparse.jar nexus_file_or_dir");
            return;
        }

        String param = args[0];

        File fileParam = new File(param);

        int goodFiles = 0;
        int totalFiles = 0;
        if(fileParam.isDirectory())
        {
           for(File file : fileParam.listFiles())
           {
                if(parseFile(file))
                    goodFiles++;

               totalFiles++;
           }
        }
        else
        {
            if(parseFile(fileParam))
                goodFiles++;

            totalFiles++;
        }

        System.out.print(String.format("%s of %s files passed", goodFiles, totalFiles));

    }

    private static boolean parseFile(File file) throws IOException {

        System.out.println(String.format("Parsing '%s'", file.getAbsolutePath()));

        try {
            Parser.parse(new FileInputStream(file));
        } catch (CoordinateParseErrorsException e)
        {
            for(CoordinateParseError error : e.Errors)
            {
               System.err.println(String.format("\t[%s, %s]: %s", error.getLineNumber(), error.getColumnNumber(), error.getMessage()));
            }
            return false;
        }

        return true;

    }
}
