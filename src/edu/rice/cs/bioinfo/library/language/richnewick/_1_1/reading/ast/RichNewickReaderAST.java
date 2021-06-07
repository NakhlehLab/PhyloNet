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

package edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.DAGFactory;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.GraphBuilder;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.csa.ASTContextAnalyser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReaderAndBuilderBase;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.csa.CSAError;
import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.io.IOException;
import java.io.InputStream;
import java.math.BigDecimal;
import java.util.Collection;
import java.util.LinkedList;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/21/13
 * Time: 5:30 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickReaderAST<G> extends RichNewickReaderAndBuilderBase<Networks, Network, GraphBuilder<G>>
{
    private final Func<RichNewickParser<Networks>> _makeParser;

    public RichNewickReaderAST(Func<RichNewickParser<Networks>> makeParser)
    {
        _makeParser = makeParser;
    }

    @Override
    protected void addNetworkToGraph(Network network, GraphBuilder<G> graphBuilder) {
        DAGFactory.makeDAG(network, graphBuilder);
    }

    @Override
    protected Collection<CSAError> detectCSAErrors(Network network, BigDecimal hybridSumTolerance) {
        return IterableHelp.toList( new ASTContextAnalyser().analyse(network, hybridSumTolerance));
    }

    @Override
    protected Iterable<Network> getNetworks(Networks networks) {
        LinkedList<Network> result = new LinkedList<Network>();

        for(NetworkNonEmpty ne : networks.Networks)
        {
            result.add(ne);
        }

        return result;
    }

    @Override
    protected Networks parse(InputStream instream) throws IOException, CoordinateParseErrorsException{
        return _makeParser.execute().parse(instream);
    }
}
