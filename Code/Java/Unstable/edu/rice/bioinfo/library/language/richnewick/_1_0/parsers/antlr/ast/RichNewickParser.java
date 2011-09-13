package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;

import edu.rice.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.bioinfo.library.language.parsing.CoordinateParseErrorDefault;
import edu.rice.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Networks;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.antlr.runtime.RecognitionException;

import java.io.*;
import java.util.LinkedList;
import java.util.List;

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

    public static Networks parse(InputStream stream) throws IOException, CoordinateParseErrorsException {
        ANTLRInputStream antlrStream = new ANTLRInputStream(stream);
        ExtendedNewickLexer lexer = new ExtendedNewickLexer(antlrStream);
        ExtendedNewickParser antlrParser = new ExtendedNewickParser(new CommonTokenStream(lexer));
        return parse(antlrParser);
    }

    static Networks parse(ExtendedNewickParser parser) throws IOException, CoordinateParseErrorsException {

        LinkedList<CoordinateParseError> errors = new LinkedList<CoordinateParseError>();
        try
        {
            parser.networks();
        }
        catch(RecognitionException e)
        {
            errors.add(new CoordinateParseErrorDefault(e.getMessage(), e.line, e.charPositionInLine));
            throw new CoordinateParseErrorsException(errors);
        }

        List<ExtendedNewickParser.ErrorWrapper> errorWraps = parser.getErrors();
        if(errors.size() > 0)
        {
            for(ExtendedNewickParser.ErrorWrapper error : errorWraps)
            {
                errors.add(new CoordinateParseErrorDefault(error.Message, error.Line, error.Col));
            }

            throw new CoordinateParseErrorsException(errors);
        }

        RuntimeException possibleException = parser.getParseStack().getException();

        if(possibleException != null)
        {
            throw possibleException;
        }

        return (Networks) parser.stack.Pop();
    }
}
