package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/31/11
 * Time: 5:44 PM
 * To change this template use File | Settings | File Templates.
 */
public interface RNewichAssignmentsBlockBody<T extends RNewickAssigmentContainer> extends Block {

    Iterable<T> getAssignments();
}
