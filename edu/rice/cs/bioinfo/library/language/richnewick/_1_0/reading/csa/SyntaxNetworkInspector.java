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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.csa;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 10:48 AM
 * To change this template use File | Settings | File Templates.
 */
public interface SyntaxNetworkInspector<N>
{
    String getNodeLabelText(N node);

    int getNodeLabelTextLineNumber(N node);

    int getNodeLabelTextColumnNumber(N node);

    String getHybridNodeIndexText(N node);

    int getHybridNodeIndexLineNumber(N node);

    int getHybridNodeIndexColumnNumber(N node);

    String getHybridNodeType(N node);

    int getHybridNodeTypeLineNumber(N node);

    int getHybridNodeTypeColumnNumber(N node);

    String getBranchLengthText(N node);

    int getBranchLengthLineNumber(N node);

    int getBranchLengthColumnNumber(N node);

    String getSupportText(N node);

    int getSupportLineNumber(N node);

    int getSupportColumnNumber(N node);

    String getProbabilityText(N node);

    int getProbabilityLineNumber(N node);

    int getProbabilityColumnNumber(N node);
}
