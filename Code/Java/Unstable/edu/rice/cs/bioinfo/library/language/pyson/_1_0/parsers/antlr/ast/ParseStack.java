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

package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.PySONNode;
import org.antlr.runtime.Token;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 1:31 PM
 * To change this template use File | Settings | File Templates.
 */
public interface ParseStack {

    public void pushIdentifier(Token ident);

    public void pushIdentList(int numElementsInList, Token startToken);

    public void pushTreesBlockBody(boolean containsTranslation);

    public void pushNetworksBlockBody(boolean containsTranslation);

    public void pushRichNewickAssignment(boolean isDefault);

    public void pushRichNewickString(String richNewickString, int line, int col);

    public void pushBlocks();

    public void pushTreeAssignment(Token treeType);

    public void pushNetworkAssignment();

    public PySONNode pop();

    public RuntimeException getException();

    public void pushPhylonetBlockBody();

    public void pushPhylonetCommandPartQuote(Token text);

    public void pushPhylonetCommandPartIdent();

    public void pushPhylonetCommandPartIdentList();

    public void pushPhylonetCommandPartSetList(Token text);

    public void pushPhylonetCommandPartIdSet(Token text);

    public void pushPhylonetCommandPartTaxaMap();

    public void pushTaxaMap(int numKeys, Token startToken);

    public void pushTaxaMapEntry(int numValues);

    public void pushPhylonetCommand(boolean includeAssigment);

    public void pushFASTAEntry(int numIdentsInDesc);

    public void pushFASTABlockBody();

    public void pushDataBlockBody(int numPairs);

}
