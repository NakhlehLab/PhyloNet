package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.Identifier;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 1:35 PM
 * To change this template use File | Settings | File Templates.
 */
public interface BlockContents {

    RichNewickAssignment getRichNewickAssigment(String identifier);

    Iterable<String> getRickNewickAssignmentIdentifiers();

    Iterable<SyntaxCommand> getCommands();


}
