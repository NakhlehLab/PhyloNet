package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;

import edu.rice.bioinfo.library.language.richnewick._1_0.ast.*;
import org.antlr.runtime.*;

import java.util.LinkedList;
import java.util.Stack;

class ParseStackAction implements ParseStack
{
    private Stack<AbstractSyntaxNode>  _parseStack = new Stack<AbstractSyntaxNode>();

    private RuntimeException _generatedException;

    public RuntimeException getException()
    {
        return _generatedException;
    }

    public void pushUnquotedText(String text)
    {
        try
        {
            _parseStack.push(new Text(text, false));
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
            _parseStack.push(new Text(insideQuotes, true));
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

    public void pushNetwork(boolean containsDescendantList)
    {
        try
        {
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

            _parseStack.push(new Network(descendants,networkInfo));
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

    public void pushBootstrap()
    {
        try
        {
            Text lengthText = (Text) _parseStack.pop();
            _parseStack.push(new BootstrapNonEmpty(lengthText));
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
                         boolean containsBootstrap, boolean containsProbability)
    {
        try
        {
            Probability probability = containsProbability ? (Probability)_parseStack.pop() : ProbabilityEmpty.Singleton;
            Bootstrap bootstrap = containsBootstrap ? (Bootstrap)_parseStack.pop() : BootstrapEmpty.Singleton;
            BranchLength branchLength = containsBranchLength ? (BranchLength)_parseStack.pop() : BranchLengthEmpty.Singleton;

            HybridNodeQualifier hybridNodeQualifier = containsHybridNodeQualifier ? (HybridNodeQualifier)_parseStack.pop() : HybridNodeQualifierEmpty.Singleton;
            NodeLabel nodeLabel = containsNodeLabel ? (NodeLabel)_parseStack.pop() : NodeLabelEmpty.Singleton;

            _parseStack.push(new NetworkInfo(nodeLabel, hybridNodeQualifier, branchLength, bootstrap, probability));
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
                text = new Text(text.Content.replace('_', ' '), text.OriginallyQuoted);
            }

            _parseStack.push(new NodeLabelNonEmpty(text));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }
    }

    public AbstractSyntaxNode Pop()
    {
        return _parseStack.pop();
    }


}
