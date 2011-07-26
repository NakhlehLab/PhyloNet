package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;

import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Network;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.antlr.runtime.RecognitionException;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 4:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickParser
{
    private RichNewickParser()
    {}

    public static Network parse(InputStream stream) throws IOException, RecognitionException
    {
        ANTLRInputStream antlrStream = new ANTLRInputStream(stream);
        ExtendedNewickLexer lexer = new ExtendedNewickLexer(antlrStream);
        ExtendedNewickParser antlrParser = new ExtendedNewickParser(new CommonTokenStream(lexer));
        return parse(antlrParser);
    }

    static Network parse(ExtendedNewickParser parser) throws IOException, RecognitionException {
        parser.network();

        RuntimeException possibleException = parser.getParseStack().getException();

        if(possibleException != null)
        {
            throw possibleException;
        }

        return (Network) parser.stack.Pop();
    }
}
