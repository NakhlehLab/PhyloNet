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

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.*;
import org.antlr.runtime.Token;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 1:31 PM
 * To change this template use File | Settings | File Templates.
 */
class ParseStackAction implements ParseStack {

    private Stack<PySONNode> _parseStack = new Stack<PySONNode>();

    private RuntimeException _exception = null;

    public RuntimeException getException() {
        return _exception;
    }

    public PySONNode pop()
    {
        return _parseStack.pop();
    }

    public void pushIdentifier(Token ident) {
        try {
            _parseStack.push(new Identifier(ident.getText(), ident.getLine(), ident.getCharPositionInLine()));
        } catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushTreeAssignment(Token treeType) {
        try {
            RichNewickAssignment assignment = (RichNewickAssignment) _parseStack.pop();

            if (PySONLexer.UTREE == treeType.getType()) {
                _parseStack.push(new TreeAssignment(false, assignment));
            } else {
                _parseStack.push(new TreeAssignment(true, assignment));
            }
        } catch (
                RuntimeException e
                )

        {
            _exception = e;
        }

    }

    public void pushNetworkAssignment() {
        try {
            RichNewickAssignment assignment = (RichNewickAssignment) _parseStack.pop();
            _parseStack.push(new NetworkAssignment(assignment));

        } catch (
                RuntimeException e
                )

        {
            _exception = e;
        }

    }

    public void pushTreesBlockBody(boolean containsTranslation) {

        try {
            if (containsTranslation) {

            } else {
                LinkedList<TreeAssignment> assignments = new LinkedList<TreeAssignment>();
                while (!_parseStack.empty() && _parseStack.peek() instanceof TreeAssignment) {
                    assignments.addFirst((TreeAssignment) _parseStack.pop());
                }

                _parseStack.push(new TreesBlockBodyWithoutTranslation(assignments));
            }

        } catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushNetworksBlockBody(boolean containsTranslation) {

        try {
            if (containsTranslation) {

            } else {
                LinkedList<NetworkAssignment> assignments = new LinkedList<NetworkAssignment>();
                while (!_parseStack.empty() && _parseStack.peek() instanceof NetworkAssignment) {
                    assignments.addFirst((NetworkAssignment) _parseStack.pop());
                }

                _parseStack.push(new NetworksBlockBody(assignments));
            }

        } catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushRichNewickAssignment(boolean isDefault) {

        try {
            RichNewickString rnString = (RichNewickString) _parseStack.pop();
            Identifier name = (Identifier) _parseStack.pop();


            _parseStack.push(new RichNewickAssignment(isDefault, name, rnString));
        } catch (RuntimeException e) {
            _exception = e;
        }

    }

    public void pushRichNewickString(String richNewickString, int line, int col) {

        try {
            _parseStack.push(new RichNewickString(richNewickString, line, col));
        } catch (RuntimeException e) {
            _exception = e;
        }
    }


    public void pushBlocks() {

        try {

            LinkedList<Block> blocks = new LinkedList<Block>();

            while (!_parseStack.isEmpty()) {
                blocks.addFirst((Block) _parseStack.pop());
            }

            _parseStack.push(new Blocks(blocks));
        } catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushPhylonetBlockBody()
    {
        try
        {
            LinkedList<PhyloNetCommand> commands = new LinkedList<PhyloNetCommand>();

            while(!_parseStack.empty() && _parseStack.peek() instanceof PhyloNetCommand)
            {
                commands.addFirst((PhyloNetCommand)_parseStack.pop());
            }

            _parseStack.push(new PhyloNetBlockBody(commands));
        }
        catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushPhylonetCommandPartQuote(Token part)
    {
        _parseStack.push(new PhyloNetCommandPartQuote(part.getText(), part.getLine(), part.getCharPositionInLine()));
    }

    public void pushPhylonetCommandPartIdent()
    {
        try
        {
            Identifier ident = (Identifier)_parseStack.pop();
            _parseStack.push(new PhyloNetCommandPartIdentifier(ident.Content, ident.Line, ident.Col));
        }
        catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushPhylonetCommandPartIdentList() {

        try
        {
            IdentList ident = (IdentList)_parseStack.pop();
            _parseStack.push(new PhyloNetCommandPartIdentList(ident, ident.Line, ident.Col));
        }
        catch (RuntimeException e) {
            _exception = e;
        }

    }

    public void pushIdentList(int numElements, Token startToken)
    {
        try
        {
            LinkedList<Identifier> elements = new LinkedList<Identifier>();
            for(int i=0; i<numElements; i++)
            {
                Identifier identifier = (Identifier)_parseStack.pop();
                elements.addAll(readIdentifier(identifier, numElements==1));
            }

            _parseStack.push(new IdentList(startToken.getLine(), startToken.getCharPositionInLine(), elements));

        }
        catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushPhylonetCommandPartSetList(Token text)
    {
        _parseStack.push(new PhyloNetCommandPartSetList(text.getText(), text.getLine(), text.getCharPositionInLine()));
    }

    public void pushPhylonetCommandPartIdSet(Token text)
    {
        _parseStack.push(new PhyloNetCommandPartIdentSet(text.getText(), text.getLine(), text.getCharPositionInLine()));
    }

    public void pushPhylonetCommandPartTaxaMap()
    {
        try
        {
            TaxaMap map = (TaxaMap)_parseStack.pop();
            _parseStack.push(new PhyloNetCommandPartTaxaMap(map));
        }
        catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushTaxaMap(int numKeys, Token startToken)
    {
        try
        {
            LinkedList<TaxaMapEntry> entries = new LinkedList<TaxaMapEntry>();

            for(int i = 0; i<numKeys; i++) {
                entries.addFirst((TaxaMapEntry)_parseStack.pop());
            }

            _parseStack.push(new TaxaMap(startToken.getLine(), startToken.getCharPositionInLine(), entries));
        }

        catch (RuntimeException e) {
            _exception = e;
        }
    }



    public void pushTaxaMapEntry(int numValues)
    {
        try
        {
            LinkedList<Identifier> values = new LinkedList<Identifier>();

            for(int i = 0; i<numValues; i++)
            {
                values.addFirst((Identifier)_parseStack.pop());

            }

            Identifier key = (Identifier)_parseStack.pop();
            _parseStack.push(new TaxaMapEntry(key, values));
        }

        catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushPhylonetCommand(boolean includeAssigment)
    {
        try
        {
            LinkedList<PhyloNetCommandPart> parts = new LinkedList<PhyloNetCommandPart>();

            while(!_parseStack.empty() && _parseStack.peek() instanceof PhyloNetCommandPart)
            {
                parts.addFirst((PhyloNetCommandPart)_parseStack.pop());
            }

            if(includeAssigment)
            {
                Identifier assignmentIdent = (Identifier) _parseStack.pop();
                _parseStack.push(new PhyloNetCommand(parts, assignmentIdent));
            }
            else
            {
                _parseStack.push(new PhyloNetCommand(parts));
            }


        }
        catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushFASTAEntry(int numIdentsInDesc)
    {
        try
        {
            Identifier sequence = (Identifier)_parseStack.pop();

            for(int i = 0; i<numIdentsInDesc; i++)
            {
                _parseStack.pop();
            }

            Identifier ident = (Identifier)_parseStack.pop();

            _parseStack.push(new FASTAEntry(ident, sequence));
        }
        catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushFASTABlockBody()
    {
        try
        {
            LinkedList<FASTAEntry> entries = new LinkedList<FASTAEntry>();
            while(!_parseStack.empty() && _parseStack.peek() instanceof FASTAEntry)
            {
                entries.add( (FASTAEntry) _parseStack.pop());
            }

            _parseStack.push(new FASTABlockBody(entries));
        }
        catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushDataBlockBody(int numPairs)
    {
        try
        {
            LinkedList<AbstractMap.SimpleImmutableEntry<Identifier,Identifier>> pairs = new LinkedList<AbstractMap.SimpleImmutableEntry<Identifier, Identifier>>();
            for(int i = 0; i<numPairs; i++)
            {
                Identifier sequence = (Identifier) _parseStack.pop();
                Identifier taxonKey = (Identifier) _parseStack.pop();
                pairs.add(new AbstractMap.SimpleImmutableEntry(taxonKey, sequence));
            }

            // nchars, ntax, and type are not currently supported.
            _parseStack.pop();
            _parseStack.pop();
            _parseStack.pop();

            _parseStack.push(new DataBlockBody(pairs));
        }
        catch (RuntimeException e) {
            _exception = e;
        }
    }


    private List<Identifier> readIdentifier(Identifier identifier, boolean checkAll){
        List<Identifier> values = new ArrayList<Identifier>();
        if(checkAll && identifier.Content.toLowerCase().equals("all")){
            for(PySONNode block: _parseStack){
                if(block instanceof TreesBlockBody){
                    TreesBlockBody tbb = (TreesBlockBody)block;
                    for(TreeAssignment ta: tbb.getAssignments()){
                        values.add(ta.getAssignment().LHSIdentifier);
                    }
                }
            }

        }

        else {

            boolean isSet = false;
            if (identifier.Content.contains("-")) {
                String[] identPairs = identifier.Content.split("-");
                char[] startingIdent = identPairs[0].toCharArray();
                int index;
                for (index = startingIdent.length - 1; index >= 0; index--) {
                    if (startingIdent[index] > '9' || startingIdent[index] < '0') {
                        break;
                    }
                }
                int startingID = Integer.parseInt(identPairs[0].substring(index + 1));
                String startingPrefix = identPairs[0].substring(0, index + 1);
                char[] endingIdent = identPairs[1].toCharArray();
                for (index = endingIdent.length - 1; index >= 0; index--) {
                    if (endingIdent[index] > '9' || endingIdent[index] < '0') {
                        break;
                    }
                }

                int endingID = Integer.parseInt(identPairs[1].substring(index + 1));
                String endingPrefix = identPairs[1].substring(0, index + 1);
                if (startingPrefix.equals(endingPrefix) && startingID <= endingID) {
                    for (int j = startingID; j <= endingID; j++) {
                        values.add(new Identifier(startingPrefix + j, identifier.Line, identifier.Col));
                    }
                    isSet = true;
                }

            }
            if (!isSet) {
                values.add(identifier);
            }
        }
        return values;
    }


}
