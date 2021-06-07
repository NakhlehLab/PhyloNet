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

package edu.rice.cs.bioinfo.library.language.hybridsimnewick._2012_6_22.printing;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.io.StringWriter;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/23/12
 * Time: 2:32 PM
 * To change this template use File | Settings | File Templates.
 */
public interface HybridSimNewickPrinter<N>
{

    void print(N root, Func1<N, Iterable<N>> getAdjacentNodes, Func1<N, Tuple<N,N>> getHybridParents, String rootBranchLength, StringWriter writer);
}