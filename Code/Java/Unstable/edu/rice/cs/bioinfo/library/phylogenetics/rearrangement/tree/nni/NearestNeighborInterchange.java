package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree.nni;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree.TreeValidator;
import edu.rice.cs.bioinfo.library.programming.Proc4;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/4/12
 * Time: 2:50 PM
 * To change this template use File | Settings | File Templates.
 */
public interface NearestNeighborInterchange<T,N,E> extends TreeValidator<N,E>
{
    void computeRearrangementsWithValidation(T tree, Proc4<T,E,E,E> rearrangementComputed);

    void computeRearrangementsWithoutValidation(T tree, Proc4<T,E,E,E> rearrangementComputed);

    public T performInterchange(T tree, E internalEdge, E swapEdgeA, E swapEdgeB);
}
