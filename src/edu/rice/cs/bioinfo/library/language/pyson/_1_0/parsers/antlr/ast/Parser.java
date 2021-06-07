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

package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorDefault;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.Blocks;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.antlr.runtime.RecognitionException;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/1/11
 * Time: 11:18 AM
 * To change this template use File | Settings | File Templates.
 */
public class Parser
{
    private Parser()
    {}

    public static Blocks parse(InputStream stream) throws IOException, CoordinateParseErrorsException {
        ANTLRInputStream antlrStream = new ANTLRInputStream(stream);
        PySONLexer lexer = new PySONLexer(antlrStream);
        PySONParser parser = new PySONParser(new CommonTokenStream(lexer));
        return parse(parser, lexer);
    }

    static Blocks parse(PySONParser parser, PySONLexer lexer) throws IOException, CoordinateParseErrorsException {

        try
        {
             parser.blocks();
        }
        catch(RecognitionException e)
        {
            CoordinateParseErrorDefault error  = new CoordinateParseErrorDefault(e.getMessage(), e.line, e.charPositionInLine);
            List<CoordinateParseError> accum = new ArrayList<CoordinateParseError>();
            accum.add(error);
            throw new CoordinateParseErrorsException(accum);
        }

        List<PySONLexer.ErrorWrapper> lexErrors = lexer.getErrors();
        if(lexErrors.size() > 0)
        {
            LinkedList<CoordinateParseError> newErrors = new LinkedList<CoordinateParseError>();

            for(PySONLexer.ErrorWrapper error : lexErrors)
            {
                newErrors.add(new CoordinateParseErrorDefault(error.Message, error.Line, error.Col));
            }

            throw new CoordinateParseErrorsException(newErrors);
        }

        List<PySONParser.ErrorWrapper> errors = parser.getErrors();
        if(errors.size() > 0)
        {
            LinkedList<CoordinateParseError> newErrors = new LinkedList<CoordinateParseError>();

            for(PySONParser.ErrorWrapper error : errors)
            {
                newErrors.add(new CoordinateParseErrorDefault(error.Message, error.Line, error.Col));
            }

            throw new CoordinateParseErrorsException(newErrors);
        }

        RuntimeException possibleException = parser.getParseStack().getException();

        if(possibleException != null)
        {
            throw possibleException;
        }

        return (Blocks) parser.stack.pop();
    }
}
