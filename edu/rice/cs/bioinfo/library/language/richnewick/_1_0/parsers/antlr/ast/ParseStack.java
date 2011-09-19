package edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.*;
import org.antlr.runtime.*;
import sun.reflect.generics.tree.VoidDescriptor;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

interface ParseStack
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


    public AbstractSyntaxNode Pop();

    public RuntimeException getException();



}
