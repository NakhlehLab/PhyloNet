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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorDefault;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickParser;
import edu.rice.cs.bioinfo.library.programming.Func;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.antlr.runtime.RecognitionException;

import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/14/11
 * Time: 4:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class ANTLRRichNewickParser implements RichNewickParser<Networks>
{
    public static final Func<RichNewickParser<Networks>> MAKE_DEFAULT_PARSER = new Func<RichNewickParser<Networks>>() {
        public ANTLRRichNewickParser execute() {
            return new ANTLRRichNewickParser();
        }
    };

    public Networks parse(InputStream stream) throws CoordinateParseErrorsException {
        try
        {
            ANTLRInputStream antlrStream = new ANTLRInputStream(stream);
            ExtendedNewickLexer lexer = new ExtendedNewickLexer(antlrStream);
            ExtendedNewickParser antlrParser = new ExtendedNewickParser(new CommonTokenStream(lexer));
            return parse(antlrParser);
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    Networks parse(ExtendedNewickParser parser) throws IOException, CoordinateParseErrorsException {

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
        if(errorWraps.size() > 0)
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

        return (Networks) parser.stack.pop();
    }
}
