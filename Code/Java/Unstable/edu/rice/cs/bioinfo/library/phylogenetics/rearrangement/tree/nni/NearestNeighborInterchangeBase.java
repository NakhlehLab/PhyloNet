package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree.nni;

import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree.TreeValidatorBase;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Proc4;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/4/12
 * Time: 4:17 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NearestNeighborInterchangeBase<T extends GraphReadOnly<N,E>, N,E> extends TreeValidatorBase<N,E> implements NearestNeighborInterchange<T,N,E>
{
    protected final Func2<N,N,E> makeEdge;

    public NearestNeighborInterchangeBase(Func2<N,N,E> makeEdge)
    {
        this.makeEdge = makeEdge;
    }

    public void computeRearrangementsWithValidation(T tree, Proc4<T,E,E,E> rearrangementComputed)
    {
        assertValidTree(tree);
        computeRearrangementsWithoutValidation(tree, rearrangementComputed);
    }
}
