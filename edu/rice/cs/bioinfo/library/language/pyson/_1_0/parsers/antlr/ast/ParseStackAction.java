package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.*;
import org.antlr.runtime.Token;
import sun.plugin.javascript.navig4.Link;

import javax.lang.model.element.NestingKind;
import java.util.LinkedList;
import java.util.Stack;

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
                    assignments.addLast((TreeAssignment) _parseStack.pop());
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
                    assignments.addLast((NetworkAssignment) _parseStack.pop());
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
                blocks.add((Block) _parseStack.pop());
            }

            _parseStack.push(new Blocks(blocks));
        } catch (RuntimeException e) {
            _exception = e;
        }
    }

    public void pushPhylonetBlockBody()
    {
        LinkedList<PhyloNetCommand> commands = new LinkedList<PhyloNetCommand>();

        while(!_parseStack.empty() && _parseStack.peek() instanceof PhyloNetCommand)
        {
            commands.addFirst((PhyloNetCommand)_parseStack.pop());
        }

        _parseStack.push(new PhyloNetBlockBody(commands));
    }

    public void pushPhylonetCommandPartQuote(Token part)
    {
        _parseStack.push(new PhyloNetCommandPartQuote(part.getText(), part.getLine(), part.getCharPositionInLine()));
    }

    public void pushPhylonetCommandPartIdent()
    {
        Identifier ident = (Identifier)_parseStack.pop();
        _parseStack.push(new PhyloNetCommandPartIdentifier(ident.Content, ident.Line, ident.Col));
    }

    public void pushPhylonetCommandPartIdentList() {

        IdentList ident = (IdentList)_parseStack.pop();
        _parseStack.push(new PhyloNetCommandPartIdentList(ident, ident.Line, ident.Col));

    }

    public void pushIdentList(int numElements, Token startToken)
    {
        LinkedList<Identifier> elements = new LinkedList<Identifier>();
        for(int i=0; i<numElements; i++)
        {
          elements.add((Identifier)_parseStack.pop());
        }

        _parseStack.push(new IdentList(startToken.getLine(), startToken.getCharPositionInLine(), elements));
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
        TaxaMap map = (TaxaMap)_parseStack.pop();
        _parseStack.push(new PhyloNetCommandPartTaxaMap(map));
    }

    public void pushTaxaMap(int numKeys, Token startToken)
    {
        LinkedList<TaxaMapEntry> entries = new LinkedList<TaxaMapEntry>();

        for(int i = 0; i<numKeys; i++)
        {
            entries.addFirst((TaxaMapEntry) _parseStack.pop());
        }

        _parseStack.push(new TaxaMap(startToken.getLine(), startToken.getCharPositionInLine(), entries));
    }

    public void pushTaxaMapEntry(int numValues)
    {
        LinkedList<Identifier> values = new LinkedList<Identifier>();

        for(int i = 0; i<numValues; i++)
        {
            values.addFirst((Identifier)_parseStack.pop());
        }

        Identifier key = (Identifier)_parseStack.pop();

        _parseStack.push(new TaxaMapEntry(key, values));
    }

    public void pushPhylonetCommand()
    {
        LinkedList<PhyloNetCommandPart> parts = new LinkedList<PhyloNetCommandPart>();

        while(!_parseStack.empty() && _parseStack.peek() instanceof PhyloNetCommandPart)
        {
            parts.addFirst((PhyloNetCommandPart)_parseStack.pop());
        }

        _parseStack.push(new PhyloNetCommand(parts));
    }

}
