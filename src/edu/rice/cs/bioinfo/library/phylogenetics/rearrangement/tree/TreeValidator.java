package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree;

import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/4/12
 * Time: 3:00 PM
 * To change this template use File | Settings | File Templates.
 */
public interface TreeValidator<N,E>
{
    void assertValidTree(GraphReadOnly<N,E> tree);
}
