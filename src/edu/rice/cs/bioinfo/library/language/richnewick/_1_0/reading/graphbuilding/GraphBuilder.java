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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;

import java.math.BigDecimal;
import java.math.BigInteger;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/26/11
 * Time: 10:52 AM
 * To change this template use File | Settings | File Templates.
 */
public interface GraphBuilder<N>
{
    public N createNode(String label);

    public N createHybridNode(String label, HybridNodeType hybridType, BigInteger hybridNodeIndex);

    public void createDirectedEdge(N tail, N tip, BigDecimal branchLength, BigDecimal support, BigDecimal probability);
}
