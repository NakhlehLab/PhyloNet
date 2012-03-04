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

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/1/11
 * Time: 2:13 PM
 * To change this template use File | Settings | File Templates.
 */

package edu.rice.cs.bioinfo.library.language.parsing;

public class CoordinateParseErrorDefault implements CoordinateParseError {

    private final String _message;

    private final int _lineNumber, _columnNumber;

    public CoordinateParseErrorDefault(String message, int lineNumber, int columnNumber)
    {
        _message = message;
        _lineNumber = lineNumber;
        _columnNumber = columnNumber;
    }

    public String getMessage() {
        return _message;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getLineNumber() {
        return _lineNumber;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getColumnNumber() {
        return _columnNumber;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
