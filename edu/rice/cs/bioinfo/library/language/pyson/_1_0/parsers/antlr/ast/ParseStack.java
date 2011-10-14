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

    public void pushPhylonetCommandPartTaxaMap(Token text);

    public void pushPhylonetCommand();

}
