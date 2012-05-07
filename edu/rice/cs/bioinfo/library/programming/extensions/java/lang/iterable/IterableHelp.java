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

package edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 5:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class IterableHelp {

    public static Object[] toArray(Iterable elements)
    {
        LinkedList<Object> accum = new LinkedList<Object>();
        for(Object element : elements)
        {
            accum.add(element);
        }

        return accum.toArray();
    }

    public static int count(Iterable elements)
    {
        int i = 0;
        for(Object element : elements)
        {
            i++;
        }

        return i;
    }
}
