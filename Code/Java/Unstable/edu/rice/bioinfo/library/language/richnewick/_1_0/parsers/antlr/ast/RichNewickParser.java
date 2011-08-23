package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;

import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReadError;
import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReadException;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Networks;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.antlr.runtime.NoViableAltException;
import org.antlr.runtime.RecognitionException;

import javax.xml.transform.ErrorListener;
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

    public static Networks parse(InputStream stream) throws IOException, RecognitionException, RichNewickReadException {
        ANTLRInputStream antlrStream = new ANTLRInputStream(stream);
        ExtendedNewickLexer lexer = new ExtendedNewickLexer(antlrStream);
        ExtendedNewickParser antlrParser = new ExtendedNewickParser(new CommonTokenStream(lexer));
        return parse(antlrParser);
    }

    static Networks parse(ExtendedNewickParser parser) throws IOException, RecognitionException, RichNewickReadException {

        parser.networks();

        List<ExtendedNewickParser.ErrorWrapper> errors = parser.getErrors();
        if(errors.size() > 0)
        {
            LinkedList<RichNewickReadError> newErrors = new LinkedList<RichNewickReadError>();

            for(ExtendedNewickParser.ErrorWrapper error : errors)
            {
                newErrors.add(new RichNewickReadError(error.Message, error.Line, error.Col));
            }

            throw new RichNewickReadException(newErrors);
        }

        RuntimeException possibleException = parser.getParseStack().getException();

        if(possibleException != null)
        {
            throw possibleException;
        }

        return (Networks) parser.stack.Pop();
    }
}
