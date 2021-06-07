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

package edu.rice.cs.bioinfo.library.language.richnewick.reading;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;

import java.io.IOException;
import java.io.InputStream;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/21/13
 * Time: 5:17 PM
 * To change this template use File | Settings | File Templates.
 */
public interface RichNewickGraphBuilder<NS,G>
{
    public RichNewickReadResult<NS> read(InputStream instream, G graphBuilder) throws CoordinateParseErrorsException, IOException;

    public RichNewickReadResult<NS> read(String richNewick, G graphBuilder) throws CoordinateParseErrorsException, IOException;

    public NS readAnyErrorToRuntimeException(InputStream instream, G graphBuilder);

    public NS readAnyErrorToRuntimeException(String richNewick, G graphBuilder);


}
