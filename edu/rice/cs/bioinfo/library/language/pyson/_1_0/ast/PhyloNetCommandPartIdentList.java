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

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 5:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloNetCommandPartIdentList extends PhyloNetCommandPart {

    public final IdentList List;

    public PhyloNetCommandPartIdentList(IdentList identList, int line, int col) {
        super(line, col);

        List = identList;
    }

    @Override
    public <R, T, E extends Exception> R execute(PhyloNetCommandPartAlgo<R, T, E> algo, T input) throws E {
        return algo.forIdenList(this, input);
    }
}
