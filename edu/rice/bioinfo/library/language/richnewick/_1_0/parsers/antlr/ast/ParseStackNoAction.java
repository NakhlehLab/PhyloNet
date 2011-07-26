package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;

import edu.rice.bioinfo.library.language.richnewick._1_0.ast.AbstractSyntaxNode;
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

    public void pushUnquotedText(String text) {
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

    public void pushNetwork(boolean containsDescendantList) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushProbability() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushBootstrap() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushBranchLength() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushNetworkInfo(boolean containsNodeLabel, boolean containsHybridNodeQualifier, boolean containsBranchLength, boolean containsBootstrap, boolean containsProbability) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void pushNodeLabel() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public AbstractSyntaxNode Pop() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
