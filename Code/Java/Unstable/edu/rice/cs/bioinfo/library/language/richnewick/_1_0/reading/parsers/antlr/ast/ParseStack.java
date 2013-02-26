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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.AbstractSyntaxNode;
import org.antlr.runtime.Token;

public interface ParseStack
{


    void pushUnquotedText(String text, int lineNumber, int columnNumber);


    void pushQuotedText(Token token);


    void pushDescendantList(int numSubtrees);


    void pushSubtree(boolean containsDescendantList);


    void pushHybridNodeQualifier(Token type, Token nodeIndex);

    void pushNetworks();

    void pushNetwork(Token rootedQualifier, boolean containsDescendantList);


    void pushProbability();


    void pushSupport();


    void pushBranchLength();


    void pushNetworkInfo(boolean containsNodeLabel, boolean containsHybridNodeQualifier, boolean containsBranchLength,
                         boolean containsSupport, boolean containsProbability);


    void pushNodeLabel();


    public AbstractSyntaxNode pop();

    public RuntimeException getException();



}
