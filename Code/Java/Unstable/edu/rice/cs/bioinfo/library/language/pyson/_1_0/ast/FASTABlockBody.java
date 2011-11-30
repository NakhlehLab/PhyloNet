package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.FASTAEntry;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.PySONNode;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/2/11
 * Time: 5:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class FASTABlockBody implements PySONNode {

    public Iterable<FASTAEntry> Entries;

    public FASTABlockBody(LinkedList<FASTAEntry> entries) {

        Entries = entries;
    }
}
