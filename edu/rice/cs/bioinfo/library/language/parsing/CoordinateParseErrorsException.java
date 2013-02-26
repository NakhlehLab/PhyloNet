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

package edu.rice.cs.bioinfo.library.language.parsing;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/1/11
 * Time: 4:45 PM
 * To change this template use File | Settings | File Templates.
 */
public class CoordinateParseErrorsException extends Exception {

    public final Iterable<CoordinateParseError> Errors;

    public CoordinateParseErrorsException(Iterable<CoordinateParseError> errors)
    {
        super(errors.iterator().next().getMessage());
        Errors = errors;
    }
}
