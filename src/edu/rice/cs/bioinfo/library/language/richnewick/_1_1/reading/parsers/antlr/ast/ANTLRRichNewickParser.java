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

package edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorDefault;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.parsers.antlr.RichNewickParserBase;
import edu.rice.cs.bioinfo.library.programming.Func;
import org.antlr.runtime.ANTLRInputStream;
import org.antlr.runtime.CommonTokenStream;
import org.antlr.runtime.RecognitionException;
import org.antlr.runtime.TokenSource;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/22/13
 * Time: 1:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class ANTLRRichNewickParser extends RichNewickParserBase<RichNewick_1_1Parser, Networks,RichNewick_1_1Parser.ErrorWrapper>
{
    public static final Func<RichNewickParser<Networks>> MAKE_DEFAULT_PARSER = new Func<RichNewickParser<Networks>>() {
           public ANTLRRichNewickParser execute() {
               return new ANTLRRichNewickParser();
           }
       };

    @Override
    protected TokenSource makeLexer(ANTLRInputStream antlrStream) {
        return new RichNewick_1_1Lexer(antlrStream);
    }

    @Override
    protected RichNewick_1_1Parser makeParser(CommonTokenStream tokenStream) {
        return new RichNewick_1_1Parser(tokenStream);
    }

    @Override
    protected Networks getGoalNonTerminal(RichNewick_1_1Parser parser) {
        return (Networks) parser.getParseStack().pop();
    }

    @Override
    protected RuntimeException getStoredRuntimeExceptionIfAny(RichNewick_1_1Parser parser) {
        return parser.getParseStack().getException();
    }

    @Override
    protected CoordinateParseError transformToCoordinateParseError(RichNewick_1_1Parser.ErrorWrapper error) {
        return new CoordinateParseErrorDefault(error.Message, error.Line, error.Col);
    }

    @Override
    protected List<RichNewick_1_1Parser.ErrorWrapper> getErrorWraps(RichNewick_1_1Parser parser) {
        return parser.getErrors();
    }

    @Override
    protected void recognizeGoalNonTerminal(RichNewick_1_1Parser parser) throws RecognitionException {
        parser.networks();
    }
}
