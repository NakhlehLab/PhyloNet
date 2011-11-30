package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.Identifier;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.PySONNode;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/2/11
 * Time: 5:24 PM
 * To change this template use File | Settings | File Templates.
 */
public class FASTAEntry implements PySONNode {

    public final Identifier Ident;

    public final Identifier Sequence;

    public FASTAEntry(Identifier ident, Identifier sequence)
    {
        Ident = ident;
        Sequence = sequence;
    }
}
