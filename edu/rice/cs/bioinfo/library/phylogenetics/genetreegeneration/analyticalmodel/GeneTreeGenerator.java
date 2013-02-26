package edu.rice.cs.bioinfo.library.phylogenetics.genetreegeneration.analyticalmodel;

import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;

import java.math.BigInteger;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/4/12
 * Time: 1:17 PM
 * To change this template use File | Settings | File Templates.
 */
public interface GeneTreeGenerator<N,E>
{
    public boolean canGenerateForNetwork(GraphReadOnly<N,E> network);

    public Iterable<GraphReadOnly> generateGeneTrees(GraphReadOnly network, BigInteger numGeneTrees);

}
