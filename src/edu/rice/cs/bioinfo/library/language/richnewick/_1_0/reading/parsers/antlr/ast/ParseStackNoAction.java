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

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/21/11
 * Time: 1:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParseStackNoAction implements ParseStack {

    public RuntimeException getException()
    {
        return null;
    }

    public void pushUnquotedText(String text, int lineNumber, int columnNumber) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushQuotedText(Token token) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushDescendantList(int numSubtrees) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushSubtree(boolean containsDescendantList) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushHybridNodeQualifier(Token type, Token nodeIndex) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushNetwork(Token rootageQualifier, boolean containsDescendantList) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushProbability() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushSupport() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushBranchLength() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushNetworkInfo(boolean containsNodeLabel, boolean containsHybridNodeQualifier, boolean containsBranchLength, boolean containsSupport, boolean containsProbability) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushNodeLabel() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushNetworks()
    {
    }

    public AbstractSyntaxNode pop() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
