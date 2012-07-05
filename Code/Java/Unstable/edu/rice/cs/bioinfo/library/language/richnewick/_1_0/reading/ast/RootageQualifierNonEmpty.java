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

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/3/11
 * Time: 5:42 PM
 * To change this template use File | Settings | File Templates.
 */
public class RootageQualifierNonEmpty implements RootageQualifier
{
    public final String Qualifier;


    public RootageQualifierNonEmpty(String qualifier)
    {
        Qualifier = qualifier;
    }

    public <R, T, E extends Exception> R execute(RootageQualifierAlgo<R, T, E> algo, T input) throws E {
        return algo.forNonEmptyQualifier(this, input);
    }

    public boolean isRooted()
    {
        return Qualifier.toLowerCase().contains("r");
    }
}
