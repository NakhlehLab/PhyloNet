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

package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 5:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class IdentList extends PySONNodeLineAndCol
{
    private final LinkedList<Identifier> _elements = new LinkedList<Identifier>();

    public final Iterable<Identifier> Elements = _elements;

    public final int ElementsCount;

    public IdentList(int line, int col, Iterable<Identifier> elements)
    {
        super(line, col);
        for(Identifier ident : elements)
        {
            _elements.add(ident);
        }

        ElementsCount = _elements.size();
    }
}
