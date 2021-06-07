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

package edu.rice.cs.bioinfo.library.language.richnewick.reading;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorDefault;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.csa.CSAError;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigDecimal;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/8/11
 * Time: 2:05 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RichNewickReaderAndBuilderBase<NS,N,G> implements RichNewickReader<NS>, RichNewickGraphBuilder<NS,G> {

    private BigDecimal _hybridSumTolerance = BigDecimal.ZERO;

    public void setHybridSumTolerance(BigDecimal tolerance)
    {
        _hybridSumTolerance = tolerance;
    }

    public RichNewickReadResult<NS> read(InputStream instream) throws CoordinateParseErrorsException, IOException {
        return read(instream, null);
    }

    public RichNewickReadResult<NS> read(String richNewick) throws CoordinateParseErrorsException, IOException {
        return read(new ByteArrayInputStream(richNewick.getBytes()));
    }

     public RichNewickReadResult<NS> read(String richNewick, G graphBuilder) throws CoordinateParseErrorsException, IOException {
        return read(new ByteArrayInputStream(richNewick.getBytes()), graphBuilder);
    }

    public NS readAnyErrorToRuntimeException(InputStream instream)
    {
        try
        {
            RichNewickReadResult<NS> result = read(instream);
            Iterator<CSAError> errors = result.getContextErrors().iterator();
            if(errors.hasNext())
            {
                throw new RuntimeException(errors.next().Message);
            }
            else
            {
                return result.getNetworks();
            }

        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }
    }

     public NS readAnyErrorToRuntimeException(String richNewick)
     {
        return readAnyErrorToRuntimeException(new ByteArrayInputStream(richNewick.getBytes()));
     }

    public NS readAnyErrorToRuntimeException(InputStream instream, G graphBuilder)
    {
       try
        {
            RichNewickReadResult<NS> result = read(instream, graphBuilder);
            Iterator<CSAError> errors = result.getContextErrors().iterator();
            if(errors.hasNext())
            {
                throw new RuntimeException(errors.next().Message);
            }
            else
            {
                return result.getNetworks();
            }

        }
        catch(Exception e)
        {
            throw new RuntimeException(e.getMessage(), e);
        }
    }

    public NS readAnyErrorToRuntimeException(String richNewick, G graphBuilder)
    {
         return readAnyErrorToRuntimeException(new ByteArrayInputStream(richNewick.getBytes()), graphBuilder);
    }

    public RichNewickReadResult<NS> read(InputStream instream, G graphBuilder) throws IOException, CoordinateParseErrorsException {

           final NS networks = parse(instream);

           final List<CSAError> errors = new LinkedList<CSAError>();

           for(N network : getNetworks(networks))
           {

               for(CSAError error : detectCSAErrors(network, _hybridSumTolerance))
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
                   addNetworkToGraph(network, graphBuilder);
               }
           }

           return new RichNewickReadResult<NS>() {
               public NS getNetworks() {
                   return networks;
               }

               public Iterable<CSAError> getContextErrors() {
                   return errors;  //To change body of implemented methods use File | Settings | File Templates.
               }
           };
    }

    protected abstract void addNetworkToGraph(N network, G graphBuilder);

    protected abstract Collection<CSAError> detectCSAErrors(N network, BigDecimal hybridSumTolerance);

    protected abstract Iterable<N> getNetworks(NS networks);

    protected abstract NS parse(InputStream instream) throws CoordinateParseErrorsException, IOException;
}
