package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;

import com.sun.org.apache.xml.internal.resolver.helpers.BootstrapResolver;
import edu.rice.bioinfo.library.language.richnewick._1_0.ast.*;
import org.antlr.runtime.*;
import sun.reflect.generics.tree.VoidDescriptor;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

interface ParseStack
{


    void pushUnquotedText(String text);


    void pushQuotedText(Token token);


    void pushDescendantList(int numSubtrees);


    void pushSubtree(boolean containsDescendantList);


    void pushHybridNodeQualifier(Token type, Token nodeIndex);


    void pushNetwork(boolean containsDescendantList);


    void pushProbability();


    void pushBootstrap();


    void pushBranchLength();


    void pushNetworkInfo(boolean containsNodeLabel, boolean containsHybridNodeQualifier, boolean containsBranchLength,
                         boolean containsBootstrap, boolean containsProbability);


    void pushNodeLabel();


    public AbstractSyntaxNode Pop();

    public RuntimeException getException();



}
