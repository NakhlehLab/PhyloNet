package edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;

import edu.rice.bioinfo.library.language.pyson._1_0.ast.*;
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

    public void pushIdentifier(String ident) {

       _parseStack.push(new Identifier(ident));
    }

    public void pushTreeAssignment(Token treeType)
    {
        RichNewickAssignment assignment = (RichNewickAssignment)_parseStack.pop();

        if(PySONLexer.UTREE == treeType.getType())
        {
            _parseStack.push(new TreeAssignment(false, assignment));
        }
        else
        {
            _parseStack.push(new TreeAssignment(true, assignment));
        }

    }

    public void pushTreesBlockBody(boolean containsTranslation) {

        if(containsTranslation)
        {

        }
        else
        {
            LinkedList<TreeAssignment> assignments = new LinkedList<TreeAssignment>();
            while(!_parseStack.empty() && _parseStack.peek() instanceof TreeAssignment)
            {
                assignments.addLast((TreeAssignment)_parseStack.pop());
            }

            _parseStack.push(new TreesBlockBodyWithoutTranslation(assignments));
        }
    }

    public void pushRichNewickAssignment(boolean isDefault) {

        Identifier name = (Identifier)_parseStack.pop();
        RichNewickString rnString = (RichNewickString)_parseStack.pop();

        _parseStack.push(new RichNewickAssignment(isDefault, name, rnString));

    }

    public void pushRichNewickString(String richNewickString) {

        _parseStack.push(new RichNewickString(richNewickString));
    }


    public void pushBlocks() {

        LinkedList<Block> blocks = new LinkedList<Block>();

        while(!_parseStack.isEmpty())
        {
            blocks.add((Block)_parseStack.pop());
        }

        _parseStack.push(new Blocks(blocks));
    }

}
