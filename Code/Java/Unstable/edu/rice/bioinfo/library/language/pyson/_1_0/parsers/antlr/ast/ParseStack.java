package edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;

import edu.rice.bioinfo.library.language.pyson._1_0.ast.PySONNode;
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

    public void pushTreesBlockBody(boolean containsTranslation);

   // public void pushNetworksBlockBody(boolean containsTranslation);

    public void pushRichNewickAssignment(boolean isDefault);

    public void pushRichNewickString(String richNewickString);

    public void pushBlocks();

    public void pushTreeAssignment(Token treeType);

    public PySONNode pop();

    public RuntimeException getException();

    public void pushPhylonetBlockBody();

    public void pushPhylonetCommandPartQuote(Token text);

    public void pushPhylonetCommandPartIdent();

    public void pushPhylonetCommand();

}
