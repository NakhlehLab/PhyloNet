package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;

import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReadError;
import edu.rice.bioinfo.library.language.richnewick._1_0.RichNewickReadException;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.Networks;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.RichNewickReaderAST;
import org.antlr.runtime.RecognitionException;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/28/11
 * Time: 2:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickReaderAST_ANTLR extends RichNewickReaderAST
{
    @Override
    protected Networks parse(InputStream instream) throws RichNewickReadException {

        try
        {
            return RichNewickParser.parse(instream);
        }
        catch (IOException e)
        {
            ArrayList<RichNewickReadError> errors = new ArrayList<RichNewickReadError>();
            errors.add(new RichNewickReadError(e.getMessage(), null, null));
            throw new RichNewickReadException(errors);
        }
        catch (RecognitionException e)
        {
            ArrayList<RichNewickReadError> errors = new ArrayList<RichNewickReadError>();
            errors.add(new RichNewickReadError(e.getMessage(), e.line, e.charPositionInLine));
            throw new RichNewickReadException(errors);
        }


    }
}
