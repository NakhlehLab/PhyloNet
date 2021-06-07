package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/4/13
 * Time: 2:55 PM
 * To change this template use File | Settings | File Templates.
 */
public interface MinSpanTreesSnpInferrer<S,E,Ex extends Exception>
{
    public Iterable<Set<E>> inferMinTrees(Set<S> sequencings) throws Ex;
}
