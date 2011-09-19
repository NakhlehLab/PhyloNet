package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.*;
import org.antlr.runtime.Token;

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

    public void pushRichNewickAssignment(boolean isDefault) {

        try {
            RichNewickString rnString = (RichNewickString) _parseStack.pop();
            Identifier name = (Identifier) _parseStack.pop();


            _parseStack.push(new RichNewickAssignment(isDefault, name, rnString));
        } catch (RuntimeException e) {
            _exception = e;
        }

    }

    public void pushRichNewickString(String richNewickString) {

        try {
            _parseStack.push(new RichNewickString(richNewickString));
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
        _parseStack.push(new PhyloNetCommandPart(part.getText(), part.getLine(), part.getCharPositionInLine()));
    }

    public void pushPhylonetCommandPartIdent()
    {
        Identifier ident = (Identifier)_parseStack.pop();
        _parseStack.push(new PhyloNetCommandPart(ident.Content, ident.Line, ident.Col));
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
