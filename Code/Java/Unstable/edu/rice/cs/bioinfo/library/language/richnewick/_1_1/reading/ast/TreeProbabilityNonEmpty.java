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

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/20/13
 * Time: 4:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class TreeProbabilityNonEmpty implements TreeProbability {

    public final String Content;

    public final String ProbString;

    public final int LineNumber, ColumnNumber;

    public TreeProbabilityNonEmpty(String content, int lineNumber, int columnNumber)
    {
        Content = content;
        LineNumber = lineNumber;
        ColumnNumber = columnNumber;

        String[] whitespaceSplit = Content.split("\\s+");
        if(whitespaceSplit[whitespaceSplit.length-1].length() == 1)
        {
            ProbString = whitespaceSplit[whitespaceSplit.length-1];
        }
        else
        {
            String part = whitespaceSplit[whitespaceSplit.length-1];
            ProbString =part.substring(0, part.length()-1);
        }

    }

    public <R, E extends Exception> R execute(TreeProbabilityAlgo<R, E> algo) {
        return algo.forNonEmpty(this);
    }
}
