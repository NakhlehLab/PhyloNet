package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/5/13
 * Time: 2:19 PM
 * To change this template use File | Settings | File Templates.
 */
public interface TransMapResult<E>
{
    public Set<Set<E>> getSolutions();
}
