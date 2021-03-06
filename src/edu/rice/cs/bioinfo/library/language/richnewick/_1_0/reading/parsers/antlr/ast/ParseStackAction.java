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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import org.antlr.runtime.Token;

import java.util.LinkedList;
import java.util.Stack;

public class ParseStackAction implements ParseStack
{
    protected Stack<AbstractSyntaxNode>  _parseStack = new Stack<AbstractSyntaxNode>();

    // often times if something is wrong with an input string the error will ultimately lead to some cast exception
    // or other runtime exception.  If we allow this runtime exception to bubble from the methods in this class
    // that error message will precede the error messages generated by ANTLR (which are more helpful).  As such
    // we can catch and store generated runtime exceptions from the methods of this class to give callers a chance
    // to instead show the ANTLR error messages first and then rethrow any exceptions from this class if desired.
    protected RuntimeException _generatedException;

    public RuntimeException getException()
    {
        return _generatedException;
    }

    public void pushUnquotedText(String text, int lineNumber, int columnNumber)
    {


        try
        {
            _parseStack.push(new Text(text, lineNumber, columnNumber, false));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }

    }

    public void pushQuotedText(Token token)
    {
        try
        {
            String text = token.getText();
            String insideQuotes = text.substring(1, text.length()-1);

            // by specification, '' char sequence INSIDE the quotes are converted to ' (ie escaping)
            insideQuotes = insideQuotes.replace("''", "'");
            _parseStack.push(new Text(insideQuotes, token.getLine(), token.getCharPositionInLine(), true));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public void pushDescendantList(int numSubtrees)
    {
        try
        {
            LinkedList<Subtree> subtress = new LinkedList<Subtree>();
            for(int i = 0; i<numSubtrees; i++)
            {
                subtress.addFirst((Subtree)_parseStack.pop());
            }

            _parseStack.push(new DescendantList(subtress));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public void pushSubtree(boolean containsDescendantList)
    {
        try
        {
            NetworkInfo info = (NetworkInfo) _parseStack.pop();

            DescendantList dl;
            if(containsDescendantList)
            {
                dl = (DescendantList) _parseStack.pop();
            }
            else
            {
                dl =  DescendantList.EMPTY_DESCENDANT_LIST;
            }



            _parseStack.push(new Subtree(dl, info));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public void pushHybridNodeQualifier(Token type, Token nodeIndex)
    {
        try
        {
            Text nodeIndexText = (Text) _parseStack.pop();
            if(type != null)
            {
                Text hybridNodeType = (Text) _parseStack.pop();
                _parseStack.push(new HybridNodeQualifierWithType(hybridNodeType, nodeIndexText));
            }
            else
            {
                _parseStack.push(new HybridNodeQualifierNonEmpty(nodeIndexText));
            }
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public void pushNetwork(Token rootageQualifier, boolean containsDescendantList)
    {
        try
        {
            RootageQualifier rq = rootageQualifier == null ? RootageQualifierEmpty.Singleton :
                                                             new RootageQualifierNonEmpty(rootageQualifier.getText());

            NetworkInfo networkInfo = (NetworkInfo)_parseStack.pop();

            DescendantList descendants;
            if(containsDescendantList)
            {
                descendants = (DescendantList) _parseStack.pop();
            }
            else
            {
                descendants = DescendantList.EMPTY_DESCENDANT_LIST;
            }

            _parseStack.push(new NetworkNonEmpty(rq, descendants,networkInfo));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public void pushNetworks()
    {
        try
        {
            LinkedList<NetworkNonEmpty> networks = new LinkedList<NetworkNonEmpty>();
            while(!_parseStack.empty())
            {
                networks.addFirst((NetworkNonEmpty) _parseStack.pop());
            }

            _parseStack.push(new Networks(networks));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public void pushProbability()
    {
        try
        {
            Text probText = (Text) _parseStack.pop();
            _parseStack.push(new ProbabilityNonEmpty(probText));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public void pushSupport()
    {
        try
        {
            Text lengthText = (Text) _parseStack.pop();
            _parseStack.push(new SupportNonEmpty(lengthText));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public void pushBranchLength()
    {
        try
        {
            Text lengthText = (Text) _parseStack.pop();
            _parseStack.push(new BranchLengthNonEmpty(lengthText));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public void pushNetworkInfo(boolean containsNodeLabel, boolean containsHybridNodeQualifier, boolean containsBranchLength,
                         boolean containsSupport, boolean containsProbability)
    {
        try
        {
            Probability probability = containsProbability ? (Probability)_parseStack.pop() : ProbabilityEmpty.Singleton;
            Support support = containsSupport ? (Support)_parseStack.pop() : SupportEmpty.Singleton;
            BranchLength branchLength = containsBranchLength ? (BranchLength)_parseStack.pop() : BranchLengthEmpty.Singleton;

            HybridNodeQualifier hybridNodeQualifier = containsHybridNodeQualifier ? (HybridNodeQualifier)_parseStack.pop() : HybridNodeQualifierEmpty.Singleton;
            NodeLabel nodeLabel = containsNodeLabel ? (NodeLabel)_parseStack.pop() : NodeLabelEmpty.Singleton;

            _parseStack.push(new NetworkInfo(nodeLabel, hybridNodeQualifier, branchLength, support, probability));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public void pushNodeLabel()
    {
        try
        {
            Text text = (Text) _parseStack.pop();

            /**
             * By specification _ chars in unquoted labels are converted to space.
             */
            if(!text.OriginallyQuoted)
            {
                text = new Text(text.Content, text.LineNumberStart, text.ColumnNumberStart, text.OriginallyQuoted);
            }

            _parseStack.push(new NodeLabelNonEmpty(text));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public AbstractSyntaxNode pop()
    {
        return _parseStack.pop();
    }


}
