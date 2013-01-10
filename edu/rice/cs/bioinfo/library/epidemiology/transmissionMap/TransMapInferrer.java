package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap;


import edu.uci.ics.jung.graph.DirectedGraph;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/19/12
 * Time: 12:42 PM
 * To change this template use File | Settings | File Templates.
 */
public interface TransMapInferrer<V,E>
{
    Set<DirectedGraph<V,E>> inferMaps();
}
