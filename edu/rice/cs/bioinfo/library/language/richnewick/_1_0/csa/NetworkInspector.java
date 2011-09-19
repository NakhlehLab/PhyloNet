package edu.rice.cs.bioinfo.library.language.richnewick._1_0.csa;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/25/11
 * Time: 3:17 PM
 * To change this template use File | Settings | File Templates.
 */
public interface NetworkInspector<N,E>
{
     Iterable<E> getAllInEdges(N node);

     N getTail(E edge);

     String getEdgeProbabilityText(E edge);
}
