package edu.rice.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 5:44 PM
 * To change this template use File | Settings | File Templates.
 */
public interface TreesBlockBody extends Block {

    Iterable<TreeAssignment> getAssignments();
}
