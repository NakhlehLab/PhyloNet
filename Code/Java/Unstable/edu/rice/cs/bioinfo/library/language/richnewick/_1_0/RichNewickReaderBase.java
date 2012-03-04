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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilder;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.graphbuilding.GraphBuilderNoAction;

import java.io.IOException;
import java.io.InputStream;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/8/11
 * Time: 2:05 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RichNewickReaderBase<T> implements RichNewickReader<T> {

    public <N> RichNewickReadResult<T> read(InputStream instream) throws CoordinateParseErrorsException, IOException {
        return read(instream, GraphBuilderNoAction.Singleton);
    }

    public abstract <N> RichNewickReadResult<T> read(InputStream instream, GraphBuilder<N> graphBuilder) throws CoordinateParseErrorsException, IOException;
}
