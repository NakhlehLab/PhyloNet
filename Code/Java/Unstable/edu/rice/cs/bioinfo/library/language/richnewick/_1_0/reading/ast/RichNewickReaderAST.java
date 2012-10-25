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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorDefault;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.csa.ASTContextAnalyser;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.csa.CSAError;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.GraphBuilder;
import edu.rice.cs.bioinfo.library.programming.Func;


import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/28/11
 * Time: 10:20 AM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickReaderAST extends RichNewickReaderBase<Networks>
{
    private final Func<RichNewickParser<Networks>> _makeParser;

    public RichNewickReaderAST(Func<RichNewickParser<Networks>> makeParser)
    {
        _makeParser = makeParser;
    }

    public RichNewickReadResult<Networks> read(InputStream instream) throws IOException, CoordinateParseErrorsException {
      return read(instream, null);
    }

    public <N> RichNewickReadResult<Networks> read(InputStream instream, GraphBuilder<N> graphBuilder) throws IOException, CoordinateParseErrorsException {

        final Networks networks = _makeParser.execute().parse(instream);

        final List<CSAError> errors = new LinkedList<CSAError>();

        for(Network network : networks.Networks)
        {

            for(CSAError error : ASTContextAnalyser.analyse(network))
            {
                errors.add(error);
            }

            if(errors.size() > 0)
            {
                LinkedList<CoordinateParseError> readErrors = new LinkedList<CoordinateParseError>();
                for(CSAError csaError : errors)
                {
                    readErrors.add(new CoordinateParseErrorDefault(csaError.Message, csaError.LineNumber, csaError.ColumnNumber) {
                    });
                }
                throw new CoordinateParseErrorsException(readErrors);
            }

            if(graphBuilder != null)
            {
                DAGFactory.makeDAG(network, graphBuilder);
            }
        }

        return new RichNewickReadResult<Networks>() {
            public Networks getNetworks() {
                return networks;
            }

            public Iterable<CSAError> getContextErrors() {
                return errors;  //To change body of implemented methods use File | Settings | File Templates.
            }
        };
    }
}
