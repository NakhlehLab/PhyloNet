/*
 * Copyright (c) 2013 Rice University.
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

package edu.rice.cs.bioinfo.library.language.richnewick.reading.parsers.antlr;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorDefault;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickParser;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.antlr.runtime.RecognitionException;
import org.antlr.runtime.TokenSource;

import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/21/13
 * Time: 3:24 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RichNewickParserBase<P,G,E> implements RichNewickParser<G>
{
    public G parse(InputStream stream) throws IOException, CoordinateParseErrorsException
    {
        ANTLRInputStream antlrStream = new ANTLRInputStream(stream);
        TokenSource lexer = makeLexer(antlrStream);

        CommonTokenStream tokenStream =  new CommonTokenStream(lexer);
        P antlrParser = makeParser(tokenStream);
        return parse(antlrParser);
    }

    protected abstract TokenSource makeLexer(ANTLRInputStream antlrStream);

    protected abstract P makeParser(CommonTokenStream tokenStream);

    private G parse(P parser) throws IOException, CoordinateParseErrorsException {

        LinkedList<CoordinateParseError> errors = new LinkedList<CoordinateParseError>();
        try
        {
            recognizeGoalNonTerminal(parser);
        }
        catch(RecognitionException e)
        {
            errors.add(new CoordinateParseErrorDefault(e.getMessage(), e.line, e.charPositionInLine));
            throw new CoordinateParseErrorsException(errors);
        }

        List<E> errorWraps = getErrorWraps(parser);
        if(errorWraps.size() > 0)
        {
            for(E error : errorWraps)
            {
                CoordinateParseError coordianteError = transformToCoordinateParseError(error);
                errors.add(coordianteError);
            }

            throw new CoordinateParseErrorsException(errors);
        }


        RuntimeException possibleException = getStoredRuntimeExceptionIfAny(parser);

        if(possibleException != null)
        {
            throw possibleException;
        }

        return getGoalNonTerminal(parser);

    }

    protected abstract G getGoalNonTerminal(P parser);

    protected abstract RuntimeException getStoredRuntimeExceptionIfAny(P parser);

    protected abstract CoordinateParseError transformToCoordinateParseError(E error);

    protected abstract List<E> getErrorWraps(P parser);

    protected abstract void recognizeGoalNonTerminal(P parser) throws RecognitionException;
}
